import csv
import os
import difflib
import statistics

import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 12
MEDIUM_SIZE = 14
LARGE_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

base_dir = "/Users/ewcss/data/ssbt/20220211_benchmark"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10"],
            "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", "MN12-L-D3(BJ)", "B97M-rV"],
            "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV"],
            "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "M06-HF", "M08-SO", "M11", "MN15", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

vac_mae = {x: dict() for x in methods}
vac_rel = {x: dict() for x in methods}
pcm_mae = {x: dict() for x in methods}
pcm_rel = {x: dict() for x in methods}

with open(os.path.join(base_dir, "abserrs_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]

        # if funct == "M06-HF":
        #     continue

        avg = float(row[-1])

        for group, functs in methods.items():
            if funct in functs:
                vac_mae[group][funct] = avg

with open(os.path.join(base_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        # if funct == "M06-HF":
        #     continue

        for group, functs in methods.items():
            if funct in functs:
                vac_rel[group][funct] = avg

# with open(os.path.join(base_dir, "abserrs_IEF-PCM.csv")) as file:
#     reader = csv.reader(file)
#     for i, row in enumerate(reader):
#         if i == 0:
#             continue
#         elif row[0].lower() == "average" or "3c" in row[0].lower():
#             continue
#         funct = row[0]
#         avg = float(row[-1])
#
#         # if funct == "M06-HF":
#         #     continue
#
#         for group, functs in methods.items():
#             if funct in functs:
#                 pcm_mae[group][funct] = avg
#
# with open(os.path.join(base_dir, "abserrs_rel_IEF-PCM.csv")) as file:
#     reader = csv.reader(file)
#     for i, row in enumerate(reader):
#         if i == 0:
#             continue
#         elif row[0].lower() == "average" or "3c" in row[0].lower():
#             continue
#         funct = row[0]
#         avg = float(row[-1])
#
#         # if funct == "M06-HF":
#         #     continue
#
#         for group, functs in methods.items():
#             if funct in functs:
#                 pcm_rel[group][funct] = avg


fig, axs = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

for i, dset in enumerate([vac_mae, vac_rel]):
    ax = axs[i]

    if i == 0:
        ax.set_ylabel("MAE (eV)")
    else:
        ax.set_ylabel("MRAE (unitless)")

    xs = ["GGA", "meta-GGA", "hybrid GGA", "hybrid meta-GGA"]
    avgs = list()
    lowlims = list()
    uplims = list()

    data = list()

    for group in xs:
        data.append(np.array(sorted(list(dset[group].values()))))

    ax.violinplot(data, [1,2,3,4], showmeans=False, showmedians=False, showextrema=False)

    quartile1 = np.zeros(4)
    medians = np.zeros(4)
    quartile3 = np.zeros(4)

    for i, d in enumerate(data):
        q1, m, q3 = np.percentile(d, [25, 50, 75])
        quartile1[i] = q1
        medians[i] = m
        quartile3[i] = q3

    whiskers = np.array([adjacent_values(sorted_array, q1, q3)
                         for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])

    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(xs)

plt.tight_layout()
fig.savefig("sp_performance_violin.png", dpi=150)

plt.show()