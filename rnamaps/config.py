"""Global plotting/style configuration.

Importing this module configures matplotlib (Agg backend, embeddable fonts)
and seaborn style. It also exposes the colour palette and line settings used
across all plots.
"""

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import seaborn as sns

sns.set_style("whitegrid", {'legend.frameon': True})

colors_dict = {
    'all': '#D3D3D3',
    'ctrl': '#408F76',
    'enh': '#F30C08',
    'sil': '#005299',
    'enhrest': '#FFB122',
    'silrest': '#6DC2F5',
    'const': '#666666',
}

linewidth = 3
dashes = False
