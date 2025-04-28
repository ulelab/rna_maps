import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import sys

# Folder containing your RNAmap TSV files
data_dir = Path("RNAmaps")
window = 600
type_of_exon = "enhanced" 

# 1) Find all files ending with "RNAmap.tsv"
files = list(data_dir.glob("*RNAmap.tsv"))
if not files:
    raise RuntimeError(f"No files ending in 'RNAmap.tsv' found in {data_dir}")

# 2) Load each file, filter for name == "enhanced", into a dict keyed by sample name
data = {}
suffix = ".MATS.JCEC_RNAmap.tsv"
for file in files:
    name = file.name
    if name.endswith(suffix):
        sample = name[:-len(suffix)]
    else:
        sample = file.stem

    df = pd.read_csv(file, sep="\t")
    df = df[df["name"] == type_of_exon]  # filter for enhanced
    if df.empty:
        continue  # skip samples with no 'enhanced' entries
    
    data[sample] = df

if not data:
    raise RuntimeError(f"No {type_of_exon} entries found in any file.")

print(f"Loaded {len(data)} samples from {len(files)} files.")

# 3) Collect and concatenate positions for each splice-site label
series_3 = [df[df.label == "middle_3ss"]["position"] for df in data.values()]
series_5 = [df[df.label == "middle_5ss"]["position"] for df in data.values()]

positions_3 = pd.concat([s for s in series_3 if not s.empty], ignore_index=True)
positions_5 = pd.concat([s for s in series_5 if not s.empty], ignore_index=True)

# Compute x-limits
min3, max3 = 0, window + 50
min5, max5 = window - 50, window * 2

# 4) Plot side-by-side
fig, (ax3, ax5) = plt.subplots(1, 2, sharey=True, figsize=(14, 6))

for sample, df in data.items():
    d3 = df[df.label == "middle_3ss"]
    if not d3.empty:
        ax3.plot(d3.position, d3["-log10pvalue_smoothed"], label=sample)
    d5 = df[df.label == "middle_5ss"]
    if not d5.empty:
        ax5.plot(d5.position, d5["-log10pvalue_smoothed"], label=sample)

# 5) Label and format
ax3.set_title(f"middle_3SS ({type_of_exon})")
ax5.set_title(f"middle_5SS ({type_of_exon})")
ax3.set_xlabel("Position")
ax5.set_xlabel("Position")
ax3.set_ylabel("-log10(pvalue_smoothed)")

ax3.set_xlim(min3, max3)
ticks3 = np.arange(min3, max3 + 1, 50)
ax3.set_xticks(ticks3)
ax3.set_xticklabels(np.arange(-1*window,  51, 50))

ax5.set_xlim(min5, max5)
ticks5 = np.arange(min5, max5 + 1, 50)
ax5.set_xticks(ticks5)
ax5.set_xticklabels(np.arange( -50, window+1, 50))
ax3.legend(loc="upper right", fontsize="small")
ax5.legend(loc="upper right", fontsize="small")

plt.tight_layout()

# 6) Save figure
fig.savefig(f"combo_rna_map_plot_{type_of_exon}.pdf")
fig.savefig(f"combo_rna_map_plot_{type_of_exon}.png", dpi=300)

plt.show()
sys.exit()
