import csv
import os
import difflib
import statistics

import numpy as np
from scipy.stats import linregress

import matplotlib.pyplot as plt

hybrid = {"global": {"PBE0": 0.25,
                     "PBE0-D3(BJ)": 0.25,
                     "B3LYP": 0.20,
                     "B3LYP-D3(BJ)": 0.20,
                     "mPW1PW91": 0.25,
                     "mPW1PW91-D3(BJ)": 0.25,
                     "M06-2X": 0.54,
                      "M06-2X-D3(0)": 0.54,
                      "M06-HF": 1.0,
                      "M06-HF-D3(0)": 1.0,
                      "M08-SO": 0.57,
                      "M08-SO-D3(0)": 0.57,
                      "MN15": 0.44,
                      "MN15-D3(0)": 0.44,
                      "BMK": 0.42,
                      "BMK-D3(BJ)": 0.42,
                      "TPSSh": 0.1,
                      "TPSSh-D3(BJ)": 0.1,
                      "SCAN0": 0.25,
                      "SCAN0-D3(BJ)": 0.25,
                      "mPWB1K": 0.44,
                      "mPWB1K-D3(BJ)": 0.44},
          "range-separated": {"LRC-wPBE": 0.0,
                              "LRC-wPBE-D3(BJ)": 0.0,
                              "LRC-wPBEh": 0.2,
                              "LRC-wPBEh-D3(BJ)": 0.2,
                              "CAM-B3LYP": 0.19,
                              "CAM-B3LYP-D3(0)": 0.19,
                              "rCAM-B3LYP": 0.18,
                              "rCAM-B3LYP-D3(0)": 0.18,
                              "HSE-HJS": 0.25,
                              "HSE-HJS-D3(BJ)": 0.25,
                              "wB97X": 0.16,
                              "wB97XD": 0.22,
                              "wB97XD3": 0.20,
                              "wB97XV": 0.17,
                              "wM06-D3": 0.27,
                              "M06-SX": 0.34,
                              "M06-SX-D3(0)": 0.34,
                              "M11": 0.43,
                              "M11-D3(0)": 0.43,
                              "revM11": 0.23,
                              "revM11-D3(0)": 0.23,
                              "wB97M-V": 0.15}
          }

SMALL_SIZE = 16
MEDIUM_SIZE = 18
LARGE_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

base_dir = "../data/results/new"

colors = {"global": "#8487B1", "range-separated": "#A9FADB"}

vac_mae = {x: dict() for x in hybrid}
vac_rel = {x: dict() for x in hybrid}
pcm_mae = {x: dict() for x in hybrid}
pcm_rel = {x: dict() for x in hybrid}

with open(os.path.join(base_dir, "abserrs_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        if "CAM" in funct:
            print(funct)

        # if funct in ["M06-HF", "M06-HF-D3(0)", "revM11-D3(0)", "TPSSh-D3(BJ)", "TPSSh"]:
        #     continue

        avg = float(row[-1])

        for group, functs in hybrid.items():
            if funct in functs:
                vac_mae[group][funct] = avg
                break

with open(os.path.join(base_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        # if funct in ["M06-HF", "M06-HF-D3(0)", "revM11-D3(0)"]:
        #     continue

        for group, functs in hybrid.items():
            if funct in functs:
                vac_rel[group][funct] = avg
                break

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
#         for group, functs in hybrid.items():
#             if funct in functs:
#                 pcm_mae[group][funct] = avg
#                 break
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
#         for group, functs in hybrid.items():
#             if funct in functs:
#                 pcm_rel[group][funct] = avg
#                 break


fig, axs = plt.subplots(2, 2, figsize=(12, 10))

for i, dset in enumerate([vac_mae, vac_rel]):
    print(i)
    bax = axs[i][0]
    ax = axs[i][1]

    if i == 0:
        bax.set_ylabel("MAE (eV)")
    else:
        bax.set_ylabel("MRAE (unitless)")

    if i == 1:
        ax.set_xlabel("HF exchange")

    gx = list()
    gy = list()
    rsx = list()
    rsy = list()

    gx_fit = list()
    gy_fit = list()

    for funct in dset["global"]:
        gx.append(hybrid["global"][funct])
        gy.append(dset["global"][funct])
        if funct not in ["M06-HF", "B3LYP"]:
            gx_fit.append(hybrid["global"][funct])
            gy_fit.append(dset["global"][funct])

    for funct in dset["range-separated"]:
        rsx.append(hybrid["range-separated"][funct])
        rsy.append(dset["range-separated"][funct])

    bp = bax.boxplot([sorted(list(dset["global"].values())),
                      sorted(list(dset["range-separated"].values()))],
                     labels=["Global", "Range-separted"],
                     patch_artist=True,
                     widths=0.6)

    for patch, color in zip(bp['boxes'], [colors["global"], colors["range-separated"]]):
        patch.set_facecolor(color)

    for median in bp['medians']:
        median.set(color='black')

    gres = linregress(np.array(gx_fit), np.array(gy_fit))
    print("Global: slope {}, intercept {}, R2 {}".format(gres.slope, gres.intercept, gres.rvalue ** 2))

    rsres = linregress(np.array(rsx), np.array(rsy))
    print("Range-separated: slope {}, intercept {}, R2 {}".format(rsres.slope, rsres.intercept, rsres.rvalue ** 2))

    fitx = np.linspace(0.1, 0.6, 100)
    fity = fitx * gres.slope + gres.intercept

    # if i == 0:
    #     ax.plot(fitx, fity, c="black")

    ax.scatter(gx, gy, c=colors["global"], edgecolors="black", alpha=0.8)
    ax.scatter(rsx, rsy, c=colors["range-separated"], edgecolors="black", alpha=0.8)

    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

plt.tight_layout(pad=3)
fig.savefig("hf_exchange_fraction.png", dpi=200)

plt.show()