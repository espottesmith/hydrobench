import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

from monty.serialization import loadfn, dumpfn

g_csvpd = loadfn('../data/g_mae_corrsvpd_opts.json')
h_csvpd = loadfn('../data/h_mae_corrsvpd_opts.json')
s_csvpd = loadfn('../data/s_mae_corrsvpd_opts.json')
g_svpd = loadfn('../data/g_mae_svpd_opts.json')
g_tzvpd = loadfn('../data/g_mae_tzvpd_opts.json')
e_tzvpd = loadfn('../data/elec_mae_tzvpd_opts.json')
method_classes = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D",
                          "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10"],
                  "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)",
                               "TPSS", "TPSS-D3(BJ)", "MN12-L",
                               "MN12-L-D3(BJ)", "B97M-rV"],
                  "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "B3LYP",
                                 "B3LYP-D3(BJ)", "CAM-B3LYP",
                                 "CAM-B3LYP-D3(0)", "mPW1PW91",
                                 "mPW1PW91-D3(BJ)", "wB97X", "wB97XD",
                                 "wB97XD3", "wB97XV"],
                  "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "M06-HF",
                                      "M08-SO", "M11", "MN15", "BMK",
                                      "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)",
                                      "SCAN0", "mPWB1K", "mPWB1K-D3(BJ)",
                                      "wB97M-V"]}
methods = list(np.concatenate([v for k, v in method_classes.items()]).flat)
sorted_g = {k: v for k, v in sorted(g_tzvpd.items(), key=lambda item: item[1])}
sorted_methods = [k for k, v in sorted_g.items()]
# x_labels = ['G_csvpd', 'H_csvpd', 'S_csvpd']
# x_labels = ['G_svpd', 'G_csvpd', 'G_tzvpd']
x_labels = ['E_tzvpd', 'G_tzvpd']
table_dict = {x: [] for x in x_labels}
for m in sorted_methods:
    for d, x in zip([e_tzvpd, g_tzvpd], x_labels): #manually change for
        # different ones
        for k, v in d.items():
            if m == k:
                table_dict[x].append(v)
info_dict = {}
for k, v in table_dict.items():
    minv = min(v)
    maxv = max(v)
    info_dict[k] = [minv, maxv]
df = pd.DataFrame(table_dict)

# A lot of the next code is to give columsn their own heatmaps,
# if we do G,H,S or something like that with different scales
df_info = pd.DataFrame(info_dict)
df_norm = pd.DataFrame()
for col in df:
    col_min = df_info[col][0]
    col_max = df_info[col][1]
    df_norm[col] = (df[col] - col_min) / (col_max - col_min)

# table = np.transpose(np.array([list(v) for k, v in table_dict.items()]))
vmin = df_norm.min().min()
vmax = df_norm.max().max()

norm_zero = (0 - vmin) / (vmax - vmin)
norm_one = (1 - vmin) / (vmax - vmin)
colors = [[0, 'darkblue'],
          [norm_zero, 'white'],
          [norm_one, 'white'],
          [1, 'darkred']
          ]
cmap = LinearSegmentedColormap.from_list('', colors, )
fig, ax = plt.subplots()
# scale the norm
df_plot = (df_norm - df_norm.min()) / (df_norm.max() - df_norm.min())

# heat map on the normalized `df_plot`
# use values in `df` to annotate
# color bar doesn't make sense so we remove it
# sns.heatmap(df_plot, annot=df, cmap='RdBu_r', cbar=False,
#             xticklabels=x_labels, yticklabels=sorted_methods)
sns.heatmap(df, annot=df, cmap='RdBu_r', cbar=False,
            xticklabels=x_labels, yticklabels=sorted_methods)
plt.show()
