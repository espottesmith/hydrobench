import csv
import os
import difflib
import statistics

import matplotlib.pyplot as plt

SMALL_SIZE = 12
MEDIUM_SIZE = 14
LARGE_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

base_dir = "../data/results/newmv"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10",],
           "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", ("MN12-L-D3(BJ)"), "B97M-rV",],
           "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "LRC-wPBE", "LRC-wPBE-D3(BJ)", "LRC-wPBEh", "LRC-wPBEh-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "rCAM-B3LYP", "rCAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "HSE-HJS", "HSE-HJS-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV",],
           "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "wM06-D3", "M06-SX", "M06-SX-D3(BJ)", "M06-HF", "M06-HF-D3(0)", "M08-SO", "M08-SO-D3(0)", "M11", "M11-D3(0)", "revM11", "revM11-D3(0)", "MN15", "MN15-D3(0)", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "SCAN0-D3(BJ)", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

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

        # if funct == "M06-HF" or funct == "M06-HF-D3(0)":
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

        # if funct == "M06-HF" or funct == "M06-HF-D3(0)" or funct == "LRC-wPBEh-D3(BJ)":            
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


fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

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
        data.append(sorted(list(dset[group].values())))
        avg = statistics.mean(dset[group].values())
        avgs.append(avg)
        group_sort = sorted(dset[group].items(), key=lambda x: x[1])
        print("\t min: {} ({}) max: {} ({}) avg: {}".format(group_sort[0][0], group_sort[0][1],
                                                    group_sort[-1][0], group_sort[-1][1],
                                                        avg))
        lowlims.append(abs(avg - group_sort[0][1]))
        uplims.append(abs(avg - group_sort[-1][1]))

    # ax.bar(x=xs, height=avgs, yerr=[lowlims, uplims], color=["#ff595e", "#ffca3a", "#8ac926", "#1982c4"])
    # ax.bar(x=range(1, len(xs) + 1), height=avgs, tick_label=xs, color=["#ff595e", "#ffca3a", "#8ac926", "#1982c4"], align="center")
    # ax.set_xticklabels(xs, rotation=30, ha="right")
    # ax2 = ax.twinx()

    bp = ax.boxplot(data, labels=xs, patch_artist=True)
    for patch, color in zip(bp['boxes'], ["#D81B60", "#FFC107", "#004D40", "#1E88E5"]):
        patch.set_facecolor(color)

    for median in bp['medians']:
        median.set(color='black')

    # ax.set_xticks(rotation=30, ha="right")
    # ax.set_xticklabels(xs, rotation=30, ha="right")

plt.tight_layout()
fig.savefig("performance_classes_mv.png", dpi=300)
# #
# plt.show()