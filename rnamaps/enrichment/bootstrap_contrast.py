"""Bootstrap-based contrast estimators for category vs. control.

For each non-control category vs. control, exons are resampled with
replacement ``B`` times and three contrasts are reported per position
with their 2.5/97.5 percentile bands:

- ``delta(p) = mean_cov_cat(p) - mean_cov_ctrl(p)`` (additive scale).
- ``log2fc(p) = log2(rate_cat(p) / rate_ctrl(p))`` (multiplicative
  scale, library-size invariant).
- ``log_odds_ratio(p) = logit(rate_cat(p)) - logit(rate_ctrl(p))``,
  where ``logit(p) = ln(p / (1 - p))``. This is the standard log
  odds-ratio contrast and is the right scale when each rate is a
  probability in ``[0, 1]`` (which is what binarised CLIP coverage
  gives). Unlike ``log2fc`` it stays bounded as ``rate -> 1`` (so
  saturated peaks don't blow up the contrast) and is symmetric under
  swapping cat <-> ctrl.

To keep ``log2fc`` and ``log_odds_ratio`` numerically stable when one
side's rate is near zero (or near one, for the logit), the rates are
regularised (see "Shrinkage" below). Three modes are supported via the
``shrinkage`` argument: ``"none"`` (default), ``"magnitude"``,
``"pseudocount"`` (legacy).

All three contrasts are computed on the per-exon coverage matrix
(binarised to 0/1 by default; see ``binarise`` in ``compute``), so each
"rate" is "fraction of exons positive at this base." That makes
``delta``, ``log2fc`` and ``log_odds_ratio`` reward categories where a
higher *proportion* of exons are covered, without rewarding raw
category size.

Shrinkage modes for ``log2fc`` / ``log_odds_ratio``
====================================================

``none`` (default)
    No regularisation. ``log2fc`` will produce ``+-inf`` at positions
    where one rate is exactly zero; ``log_odds_ratio`` clips rates to
    ``[eps_safe, 1 - eps_safe]`` only so that ``logit`` stays finite
    (otherwise the same: no shrinkage). The default is ``none``
    because the previous magnitude default with ``tau = 0.05`` was
    calibrated for binarised CLIP data only -- it makes no sense for
    continuous inputs (e.g. AI prediction scores via
    ``binarise=False``) where the "rate" axis is not in ``[0, 1]``.
    Switch to ``magnitude`` explicitly if you want shrinkage at sparse
    positions.

``magnitude``
    The log2 fold change (and the log odds ratio) are pulled toward
    zero by a weight that depends on the *larger* of the two rates --
    not on the number of exons. Specifically

        raw     = log2((cat + eps_safe) / (ctrl + eps_safe))
        scale   = max(cat, ctrl)
        weight  = tau / (scale + tau)
        shrunk  = (1 - weight) * raw

    Using the larger rate is the key design choice. It says "if at
    least one side has a substantial rate (above ``tau``), this is a
    real signal -- don't shrink. If *both* sides are tiny, the
    apparent fold change is unreliable -- shrink toward zero."

    - Both rates small in absolute terms (``max < tau``):
      ``weight -> 1`` and ``shrunk -> 0`` (heavy shrinkage). This is
      the case the user usually wants tamed: a single-exon flip
      against a near-zero baseline.
    - At least one rate well above ``tau``: ``weight -> 0`` and the
      contrast is essentially the raw log2 fold change. Real peaks
      (high cat rate vs low ctrl, or vice versa) are preserved.
    - One side exactly zero: still bounded (``eps_safe`` inside the
      log keeps ``raw`` finite); the magnitude weight then decides
      whether to keep it or shrink it based on the other side.

    The transform depends only on the rate magnitudes at each
    position -- not on the number of category or control exons.
    ``tau`` defaults to ``0.05`` (5% rate), a sensible "this is a
    small rate" threshold for binarised CLIP coverage. ``eps_safe``
    defaults to ``1e-3``, large enough to keep ``log2`` numerically
    well-behaved at exact zeros and small enough not to bias any
    rate well above background. Override ``tau`` with the
    ``shrinkage_scale`` argument.

``pseudocount`` (legacy)
    Adds a small additive epsilon to both rates before the log ratio.
    Default ``eps = max(1e-3, pseudocount_frac * median(rate_ctrl))``;
    override with a fixed ``pseudocount`` value. For
    ``log_odds_ratio`` the same ``eps`` is used to clip rates into
    ``[eps, 1 - eps]`` before the logit.

The bootstrap reflects sampling uncertainty in the estimator. It does
*not* correct for bias from contaminated controls (silently regulated
exons in the control pool). Use the ``--control_set`` flag for that.

When ``smoothing > 1``, each bootstrap iteration's per-position contrast
is convolved with a centred Gaussian-weighted rolling mean (matching the
kernel in ``rnamaps.permutation._smooth``) *before* taking the across-
iteration mean and percentiles. Smoothing first makes the CI band the
correct uncertainty band for the smoothed estimator; smoothing the
percentiles afterwards would understate uncertainty at sharp features.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from rnamaps.coverage import aggregate_legacy_columns
from rnamaps.enrichment._base import EnrichmentResult
from rnamaps.permutation import _build_coverage_matrix


# Allowed values for the ``shrinkage`` argument of :func:`compute`.
SHRINKAGE_MODES = ("magnitude", "pseudocount", "none")

# Default shrinkage scale (rate threshold below which heavy shrinkage
# kicks in). 0.05 = 5% rate. Chosen so that a position must have at
# least one side covered by >~5% of exons to count as a "real" rate
# worth not shrinking; below that, log2fc is pulled toward zero.
# Override with the ``shrinkage_scale`` argument.
_DEFAULT_SHRINKAGE_SCALE = 0.05

# Additive constant inside the log so the ratio is finite when one rate
# is exactly zero. 1e-3 is much smaller than any biologically plausible
# rate worth preserving (so it doesn't bias real peaks) but large
# enough to keep ``log2`` well-behaved at zero-rate positions.
_LOG_SAFETY_EPSILON = 1e-3


def _smooth_along_positions(
    arr: np.ndarray, smoothing: int
) -> np.ndarray:
    """Centred Gaussian rolling-mean smoothing along the position axis.

    Matches ``rnamaps.permutation._smooth`` (window ``smoothing``, Gaussian
    weights with ``std=2``). Accepts a 1D ``(n_pos,)`` array or a 2D
    ``(n_boot, n_pos)`` array; for the 2D case each bootstrap row is
    smoothed independently in vectorised C via ``DataFrame.rolling``.

    Edge positions where the window does not fully fit become NaN, mirroring
    the existing permutation-test behaviour.
    """
    if smoothing is None or smoothing <= 1:
        return arr
    if arr.ndim == 1:
        s = pd.Series(arr)
        return s.rolling(
            smoothing, center=True, win_type='gaussian'
        ).mean(std=2).to_numpy()
    # 2D: transpose so positions are rows and each bootstrap iteration is a
    # column; pandas rolls along axis 0 (the default) per-column → smooths
    # each bootstrap row in vectorised C. Transpose back at the end.
    df = pd.DataFrame(arr.T)
    smoothed = df.rolling(
        smoothing, center=True, win_type='gaussian'
    ).mean(std=2).to_numpy()
    return smoothed.T


def _bootstrap_means(
    matrix: np.ndarray,
    n: int,
    n_boot: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Bootstrap mean coverage per position from `matrix`.

    Parameters
    ----------
    matrix : ndarray, shape (n_rows, n_positions)
        Per-exon coverage rows.
    n : int
        Number of rows to draw with replacement per iteration (usually
        ``matrix.shape[0]``).
    n_boot : int
    rng : np.random.Generator

    Returns
    -------
    ndarray, shape (n_boot, n_positions)
    """
    if matrix.ndim != 2 or matrix.shape[0] == 0:
        n_pos = matrix.shape[1] if matrix.ndim == 2 else 0
        return np.zeros((n_boot, n_pos), dtype=np.float64)
    out = np.empty((n_boot, matrix.shape[1]), dtype=np.float64)
    n_rows = matrix.shape[0]
    for b in range(n_boot):
        idx = rng.integers(0, n_rows, size=n)
        out[b] = matrix[idx].mean(axis=0)
    return out


def _safe_logit(p: np.ndarray, eps_safe: float = _LOG_SAFETY_EPSILON) -> np.ndarray:
    """``ln(p / (1 - p))`` with ``p`` clipped into ``[eps, 1 - eps]``.

    The clip keeps the logit finite for ``p = 0`` and ``p = 1`` exactly
    (which is common with binarised coverage when zero or all of the
    exons in a category have signal at a position) without changing
    values away from the boundary. ``eps_safe = 1e-3`` corresponds to
    a maximal |logit| of about ``ln(999) ~ 6.9``, well above any
    biologically interesting odds.
    """
    safe = np.clip(p, eps_safe, 1.0 - eps_safe)
    return np.log(safe / (1.0 - safe))


def _shrink_log_odds_by_magnitude(
    cat_rate: np.ndarray,
    ctrl_rate: np.ndarray,
    tau: float,
    eps_safe: float = _LOG_SAFETY_EPSILON,
) -> np.ndarray:
    """Shrink ``logit(cat) - logit(ctrl)`` by the same weight used for
    log2fc magnitude shrinkage.

    ``weight = tau / (max(cat, ctrl) + tau)`` -- exactly mirrors
    :func:`_shrink_log2fc_by_magnitude` so that both contrasts shrink
    at the same positions for the same reason (both rates tiny =
    unreliable contrast, pulled toward zero).
    """
    cat_arr = np.asarray(cat_rate, dtype=np.float64)
    ctrl_arr = np.asarray(ctrl_rate, dtype=np.float64)
    raw = _safe_logit(cat_arr, eps_safe) - _safe_logit(ctrl_arr, eps_safe)
    scale = np.maximum(cat_arr, ctrl_arr)
    weight = tau / (scale + tau)
    return (1.0 - weight) * raw


def _shrink_log2fc_by_magnitude(
    cat_rate: np.ndarray,
    ctrl_rate: np.ndarray,
    tau: float,
    eps_safe: float = _LOG_SAFETY_EPSILON,
) -> np.ndarray:
    """Shrink ``log2(cat / ctrl)`` toward 0 by a weight that depends on
    the *larger* of the two rates.

    Concretely

        raw     = log2((cat + eps_safe) / (ctrl + eps_safe))
        scale   = max(cat, ctrl)
        weight  = tau / (scale + tau)
        shrunk  = (1 - weight) * raw

    Using ``max(cat, ctrl)`` is the key choice: if at least one side
    has a substantial rate (well above ``tau``) we have real signal,
    so the contrast is preserved. If *both* sides are tiny, the
    apparent fold change is dominated by noise / a single-exon flip
    and the contrast is pulled toward zero.

    Properties:

    * ``max < tau`` (both rates tiny) -> ``weight -> 1`` ->
      ``shrunk -> 0`` (heavy shrinkage).
    * ``max >> tau`` (at least one side has real signal) ->
      ``weight -> 0`` -> ``shrunk -> raw`` (no shrinkage). This
      preserves real peaks even when the *other* side happens to be
      zero, which the geometric-mean variant of this formula does not.
    * Equal rates -> ``raw = 0`` -> ``shrunk = 0`` (trivially).
    * One side exactly zero -> ``raw`` is bounded by ``log2(cat / eps_safe)``
      so no ``+-inf`` spikes. The weight then decides: if the other
      side is big, signal is preserved; if both are small, shrinkage
      kicks in.
    * **Independent of the number of exons** -- the decision depends
      only on the rate values at this position.

    ``eps_safe`` keeps ``log2`` finite at exact zeros. It is small
    enough not to bias any biologically meaningful rate downward.
    """
    cat_arr = np.asarray(cat_rate, dtype=np.float64)
    ctrl_arr = np.asarray(ctrl_rate, dtype=np.float64)
    safe_cat = cat_arr + eps_safe
    safe_ctrl = ctrl_arr + eps_safe
    raw = np.log2(safe_cat / safe_ctrl)
    scale = np.maximum(cat_arr, ctrl_arr)
    weight = tau / (scale + tau)
    return (1.0 - weight) * raw


def compute(
    region,
    exon_categories: pd.Series,
    label: str,
    *,
    rng: np.random.Generator,
    n_boot: int = 1000,
    shrinkage: str = "none",
    shrinkage_scale: Optional[float] = None,
    pseudocount: Optional[float] = None,
    pseudocount_frac: float = 0.01,
    bootstrap_control_fixed: bool = False,
    smoothing: int = 1,
    control_label: str = "control",
    ci_low: float = 2.5,
    ci_high: float = 97.5,
    binarise: bool = True,
) -> EnrichmentResult:
    """Compute bootstrap-CI contrasts (delta and log2fc) for each non-control category.

    Parameters
    ----------
    region : rnamaps.coverage.RegionCoverage or pd.DataFrame
        Per-exon coverage for one splice-site region. Long-form
        DataFrame input is still accepted for backward compat with
        existing tests.
    exon_categories : Series
        Counts per category.
    label : str
        Splice-site region label, e.g. ``"middle_3ss"``.
    rng : np.random.Generator
    n_boot : int
        Bootstrap iterations.
    shrinkage : {"none", "magnitude", "pseudocount"}
        Regularisation applied to ``log2fc`` and ``log_odds_ratio``.
        ``"none"`` (default) does no shrinkage: ``log2fc`` will
        produce ``+-inf`` at zero-rate positions, ``log_odds_ratio``
        clips rates into ``[eps_safe, 1 - eps_safe]`` only to keep
        the logit finite (no other shrinkage). ``"magnitude"`` shrinks
        both contrasts toward zero by a weight that depends on the
        *larger* of the two rates -- not on the number of exons. When
        both rates are small in absolute terms, the contrasts are
        pulled to zero; when at least one rate is well above the
        shrinkage scale, the raw contrast is preserved.
        ``"pseudocount"`` adds an additive epsilon to both rates
        (legacy behaviour) and, for ``log_odds_ratio``, clips into
        ``[eps, 1 - eps]`` before the logit.
    shrinkage_scale : float or None
        Used only with ``shrinkage="magnitude"``. Sets ``tau`` -- the
        rate scale at which shrinkage transitions from heavy to light.
        ``None`` (default) means ``tau = 0.01`` (1% rate), a sensible
        "this is a small rate" threshold for binarised CLIP coverage.
        Pass a smaller value for gentler shrinkage that only bites at
        very small rates; pass a larger value for more aggressive
        shrinkage even at moderate rates.
    pseudocount : float or None
        Used only when ``shrinkage="pseudocount"``. If not None, fixed
        pseudocount for log2fc. Otherwise an adaptive
        ``eps = max(1e-3, pseudocount_frac * median(rate_ctrl over region))``
        is used per region/category.
    pseudocount_frac : float
        Fraction of regional median ctrl rate when adaptive pseudocount
        is in use.
    bootstrap_control_fixed : bool
        If True, treat the control mean as a constant and skip resampling
        the control group. Equivalent up to negligible variance when
        ``n_ctrl >= 2000`` and ~5-10x faster.
    smoothing : int
        Centred Gaussian rolling-mean window (in positions) applied to
        each bootstrap iteration's ``delta`` / ``log2fc`` before taking
        the across-iteration mean and percentiles. ``smoothing <= 1``
        disables smoothing. Matches the kernel used by ``permutation_z``.
    control_label : str
    ci_low, ci_high : float
        Percentile bounds for the CI band.
    binarise : bool, default True
        Whether to threshold the per-exon coverage matrix to 0/1 before
        bootstrapping. The default matches CLIP semantics where each
        "rate" is the fraction of exons positive at that base, which is
        what makes ``delta``, ``log2fc`` and the magnitude-shrinkage
        scale comparable across positions. Set to ``False`` for
        continuous inputs (e.g. AI prediction scores or density
        tracks); ``delta`` then becomes a mean-signal difference and
        ``log2fc`` a ratio of mean signal levels. Note that the
        magnitude-shrinkage default ``tau = 0.05`` assumes a [0, 1]
        rate axis -- pass ``shrinkage_scale`` matched to the typical
        scale of your continuous signal when ``binarise=False``.

    Returns
    -------
    EnrichmentResult
        ``plot_df`` carries ``delta``, ``delta_lo``, ``delta_hi``,
        ``log2fc``, ``log2fc_lo``, ``log2fc_hi``, ``log_odds_ratio``,
        ``log_odds_ratio_lo``, ``log_odds_ratio_hi`` plus the legacy
        coverage / fold_change parity columns. Shrinkage metadata is
        attached as the ``shrinkage``, ``tau`` and ``pseudocount``
        columns (the unused ones are filled with NaN for the chosen
        mode).
    """
    if shrinkage not in SHRINKAGE_MODES:
        raise ValueError(
            f"shrinkage must be one of {SHRINKAGE_MODES}, got {shrinkage!r}"
        )
    if control_label not in exon_categories.index:
        raise ValueError(
            f"bootstrap_contrast requires a '{control_label}' category."
        )

    categories = [c for c in exon_categories.index if c != control_label]
    n_ctrl_total = int(exon_categories.loc[control_label])

    agg = aggregate_legacy_columns(region, exon_categories, control_label)

    plot_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            region, cat, control_label, binarise=binarise,
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl_total == 0:
            logging.warning(
                f"[bootstrap_contrast] Skipping {cat} ({label}): empty matrix."
            )
            continue

        cat_mat = matrix[is_cat]
        ctrl_mat = matrix[~is_cat]
        n_c = cat_mat.shape[0]
        n_ctrl = ctrl_mat.shape[0]

        if not bootstrap_control_fixed and n_c > 0 and n_ctrl >= 20 * n_c:
            logging.info(
                f"[bootstrap_contrast] {cat} ({label}): n_ctrl={n_ctrl} "
                f">= 20x n_{cat}={n_c}. --bootstrap_control_fixed would be "
                f"~{n_ctrl / n_c:.0f}x faster with negligible CI change."
            )

        ctrl_mean = ctrl_mat.mean(axis=0) if n_ctrl > 0 else np.zeros(
            matrix.shape[1]
        )

        # Resolve shrinkage hyper-parameters once per (region, category).
        eps_ = np.nan
        tau_ = np.nan
        if shrinkage == "pseudocount":
            if pseudocount is not None:
                eps_ = float(pseudocount)
            else:
                med = float(np.median(ctrl_mean)) if n_ctrl > 0 else 0.0
                eps_ = max(1e-3, pseudocount_frac * med)
        elif shrinkage == "magnitude":
            tau_ = (float(shrinkage_scale)
                    if shrinkage_scale is not None
                    else _DEFAULT_SHRINKAGE_SCALE)
            logging.info(
                f"[bootstrap_contrast] {cat} ({label}): "
                f"magnitude shrinkage tau={tau_:.4g} (log2fc shrunk "
                f"toward 0 where geo-mean of rates << tau)."
            )

        cat_means = _bootstrap_means(cat_mat, n_c, n_boot, rng)
        if bootstrap_control_fixed:
            ctrl_means = np.broadcast_to(
                ctrl_mean, (n_boot, ctrl_mean.size)
            )
        else:
            ctrl_means = _bootstrap_means(ctrl_mat, n_ctrl, n_boot, rng)

        delta_b = cat_means - ctrl_means

        if shrinkage == "magnitude":
            log2fc_b = _shrink_log2fc_by_magnitude(
                cat_means, ctrl_means, tau_,
            )
            log_odds_b = _shrink_log_odds_by_magnitude(
                cat_means, ctrl_means, tau_,
            )
        elif shrinkage == "pseudocount":
            log2fc_b = np.log2((cat_means + eps_) / (ctrl_means + eps_))
            # Clip rates into [eps_, 1 - eps_] before the logit so the
            # additive pseudocount has a comparable effect to the
            # log-ratio case (bounding the contrast away from infinity
            # by the same magnitude).
            log_odds_b = (_safe_logit(cat_means, eps_safe=eps_)
                          - _safe_logit(ctrl_means, eps_safe=eps_))
        else:  # "none"
            with np.errstate(divide='ignore', invalid='ignore'):
                log2fc_b = np.log2(cat_means / ctrl_means)
            # Even in "none" mode we still apply the minimal logit
            # safety clip; otherwise rates of exactly 0 or 1 yield
            # ``-inf`` / ``+inf`` everywhere they occur and the
            # percentile bands degenerate. The default eps (1e-3)
            # caps |log_odds_ratio| at ~6.9, well above any
            # biologically meaningful effect, so this clip is
            # essentially invisible at non-saturating rates.
            log_odds_b = (_safe_logit(cat_means)
                          - _safe_logit(ctrl_means))

        # Smooth each bootstrap iteration along the position axis before
        # collapsing to mean/percentiles so the CI band is the uncertainty
        # of the smoothed estimator, not the smoothed uncertainty.
        delta_b = _smooth_along_positions(delta_b, smoothing)
        log2fc_b = _smooth_along_positions(log2fc_b, smoothing)
        log_odds_b = _smooth_along_positions(log_odds_b, smoothing)

        cat_df = pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            'delta': np.nanmean(delta_b, axis=0),
            'delta_lo': np.nanpercentile(delta_b, ci_low, axis=0),
            'delta_hi': np.nanpercentile(delta_b, ci_high, axis=0),
            'log2fc': np.nanmean(log2fc_b, axis=0),
            'log2fc_lo': np.nanpercentile(log2fc_b, ci_low, axis=0),
            'log2fc_hi': np.nanpercentile(log2fc_b, ci_high, axis=0),
            'log_odds_ratio': np.nanmean(log_odds_b, axis=0),
            'log_odds_ratio_lo': np.nanpercentile(log_odds_b, ci_low, axis=0),
            'log_odds_ratio_hi': np.nanpercentile(log_odds_b, ci_high, axis=0),
            'shrinkage': shrinkage,
            'tau': tau_,
            'pseudocount': eps_,
        })
        plot_rows.append(cat_df)

    ctrl_positions = agg.loc[agg['name'] == control_label, 'position'].to_numpy()
    if len(ctrl_positions):
        plot_rows.append(pd.DataFrame({
            'name': control_label,
            'position': ctrl_positions,
            'label': label,
            'delta': 0.0,
            'delta_lo': 0.0,
            'delta_hi': 0.0,
            'log2fc': 0.0,
            'log2fc_lo': 0.0,
            'log2fc_hi': 0.0,
            'log_odds_ratio': 0.0,
            'log_odds_ratio_lo': 0.0,
            'log_odds_ratio_hi': 0.0,
            'shrinkage': shrinkage,
            'tau': np.nan,
            'pseudocount': np.nan,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label',
            'delta', 'delta_lo', 'delta_hi',
            'log2fc', 'log2fc_lo', 'log2fc_hi',
            'log_odds_ratio', 'log_odds_ratio_lo', 'log_odds_ratio_hi',
            'shrinkage', 'tau', 'pseudocount',
        ])
    else:
        plot_df = pd.concat(plot_rows, ignore_index=True)

    plot_df = plot_df.merge(
        agg[['name', 'position', 'coverage', 'number_exons', 'norm_coverage',
             'control_coverage', 'control_number_exons',
             'control_norm_coverage', 'fold_change']],
        on=['name', 'position'], how='left'
    )

    return EnrichmentResult(
        plot_df=plot_df,
        clusters_df=pd.DataFrame(),
        plot_kind='ribbon',
        y_columns=('delta', 'log2fc', 'log_odds_ratio'),
        method_name='bootstrap_contrast',
        ylabel='bootstrap contrast vs control',
    )
