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

base_dir = "../data/results/new"

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
    "CAM-B3LYP": {"D3": "CAM-B3LYP-D3(BJ)"},
    "rCAM-B3LYP": {"D3": "rCAM-B3LYP-D3(BJ)"},
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
    "M06-SX": {"D3": "M06-SX-D3(0)"},
    "M11": {"D3": "M11-D3(0)"},
    "revM11": {"D3": "revM11-D3(0)"},
}

colors = {"D2": "#23022E", "D3": "#573280", "VV10": "#F896D8"}

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


fig, axs = plt.subplots(1, 2, figsize=(8, 8))

for i, dset in enumerate([vac_mae, vac_rel]):
    ax = axs[i]

    if i % 2 == 0:
        ax.set_xlabel("Base MAE (eV)")
    else:
        ax.set_xlabel("Base MRAE (unitless)")

    ax.set_ylabel("Relative error (unitless)")

    d2_xs = list()
    d3_xs = list()
    vv_xs = list()

    d2_ys = list()
    d3_ys = list()
    vv_ys = list()

    ax.set_xticks(xs, x_labels, rotation=30, ha="right")

plt.tight_layout()
# plt.legend()

fig.savefig("sp_dispersion_effect.png", dpi=150)
# plt.show()