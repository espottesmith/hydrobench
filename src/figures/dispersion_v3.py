import csv
import os
import difflib
import statistics

from collections import defaultdict

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

base_dir = "../data/results/new"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10",],
           "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", ("MN12-L-D3(BJ)"), "B97M-rV",],
           "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "LRC-wPBE", "LRC-wPBE-D3(BJ)", "LRC-wPBEh", "LRC-wPBEh-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "rCAM-B3LYP", "rCAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "HSE-HJS", "HSE-HJS-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV",],
           "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "wM06-D3", "M06-SX", "M06-SX-D3(BJ)", "M06-HF", "M06-HF-D3(0)", "M08-SO", "M08-SO-D3(0)", "M11", "M11-D3(0)", "revM11", "revM11-D3(0)", "MN15", "MN15-D3(0)", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "SCAN0-D3(BJ)", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

re_bases = {
    "PBE": {"D3": "PBE-D3(BJ)"},
    "BLYP": {"D3": "BLYP-D3(BJ)"},
    "mPW91": {"D3": "mPW91-D3(BJ)"},
    "M06-L": {"D3": "M06-L-D3(0)"},
    "SCAN": {"D3": "SCAN-D3(BJ)"},
    "TPSS": {"D3": "TPSS-D3(BJ)"},
    "MN12-L": {"D3": "MN12-L-D3(BJ)"},
    "PBE0": {"D3": "PBE0-D3(BJ)"},
    "B3LYP": {"D3": "B3LYP-D3(BJ)"},
    "mPW1PW91": {"D3": "mPW1PW91-D3(BJ)"},
    "LRC-wPBE": {"D3": "LRC-wPBE-D3(BJ)"},
    "LRC-wPBEh": {"D3": "LRC-wPBEh-D3(BJ)"},
    "CAM-B3LYP": {"D3": "CAM-B3LYP-D3(0)"},
    "rCAM-B3LYP": {"D3": "rCAM-B3LYP-D3(0)"},
    "HSE-HJS": {"D3": "HSE-HJS-D3(BJ)"},
    "wB97X": {"D2": "wB97XD", "D3": "wB97XD3", "VV10": "wB97XV"},
    "M06-2X": {"D3": "M06-2X-D3(0)"},
    "M06-HF": {"D3": "M06-HF-D3(0)"},
    "M08-SO": {"D3": "M08-SO-D3(0)"},
    "MN15": {"D3": "MN15-D3(0)"},
    "BMK": {"D3": "BMK-D3(BJ)"},
    "TPSSh": {"D3": "TPSSh-D3(BJ)"},
    "SCAN0": {"D3": "SCAN0-D3(BJ)"},
    "mPWB1K": {"D3": "mPWB1K-D3(BJ)"},
    "M06-SX": {"D3": "M06-SX-D3(BJ)"},
    "M11": {"D3": "M11-D3(0)"},
    "revM11": {"D3": "revM11-D3(0)"},
}

shapes = {"D2": "v", "D3": "o", "VV10": "P"}
colors = {"GGA": "#D81B60", "meta-GGA": "#FFC107", "hybrid GGA": "#004D40", "hybrid meta-GGA": "#1E88E5"}

lims = [(0.03, 0.2), (0.1, 0.8)]

vac_mae = dict()
vac_rel = dict()
pcm_mae = dict()
pcm_rel = dict()

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
        vac_mae[funct] = avg

with open(os.path.join(base_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]

        # if funct == "M06-HF" or funct == "M06-HF-D3(0)" or funct == "LRC-wPBEh-D3(BJ)":            
        #     continue                                         

        avg = float(row[-1])
        vac_rel[funct] = avg


fig, axs = plt.subplots(2, 1, figsize=(8, 8))

for i, dset in enumerate([vac_mae, vac_rel]):
    ax = axs[i]

    if i % 2 == 0:
        ax.set_xlabel("Base MAE (eV)")
    else:
        ax.set_xlabel("Base MRAE (unitless)")

    ax.set_ylabel("Change in error (eV)")

    xs = {"D2": defaultdict(list), "D3": defaultdict(list), "VV10": defaultdict(list)}
    ys = {"D2": defaultdict(list), "D3": defaultdict(list), "VV10": defaultdict(list)}

    for base, subs in re_bases.items():
        if i == 0:
            this_data = vac_mae
        else:
            this_data = vac_rel

        this_group = None
        for gname, group in methods.items():
            if base in group:
                this_group = gname

        base_error = this_data.get(base)
        for disp in ["D2", "D3", "VV10"]:
            if disp in subs:
                disp_error = this_data.get(subs[disp])
                try:                    
                    diff_error = disp_error - base_error
                    xs[disp][this_group].append(base_error)
                    ys[disp][this_group].append(diff_error)
                    print(base, disp, group, diff_error)
                except TypeError:
                    print(base, subs[disp])

    print()
    ax.hlines(0.0, lims[i][0], lims[i][1], colors='k', linestyle='dashed')
    ax.set_xlim(lims[i][0], lims[i][1])

    if i == 0:
        ax.set_xticks([0.05, 0.1, 0.15, 0.2])

    for group in methods:
        ax.scatter(xs["D2"][group], ys["D2"][group], c=colors[group], edgecolors="black", alpha=0.8, label=f"{group}-D2", s=60, marker=shapes["D2"])
        ax.scatter(xs["D3"][group], ys["D3"][group], c=colors[group], edgecolors="black", alpha=0.8, label=f"{group}-D3", s=60, marker=shapes["D3"])
        ax.scatter(xs["VV10"][group], ys["VV10"][group], c=colors[group], edgecolors="black", alpha=0.8, label=f"{group}-VV10", s=60, marker=shapes["VV10"])

plt.tight_layout()
# plt.legend()

fig.savefig("sp_dispersion_effect_v3.png", dpi=200)
# plt.show()