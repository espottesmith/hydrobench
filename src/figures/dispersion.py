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

base_dir = "/Users/ewcss/data/ssbt/20220211_benchmark"

exclude = ["VV10", "rVV10", "B97M-rV",  "M06-HF", "M08-SO", "M11", "MN15", "SCAN0", "wB97M-V"]

colors = {"None": "#f72585", "D2": "#b5179e", "D3": "#7209b7", "VV10": "#3a0ca3"}

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

        if funct in exclude:
            continue

        else:
            if funct.endswith("D"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D"):
                    base = base[:-1]
                disp = "D2"
            elif funct.endswith("D3"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D3"):
                    base = base[:-2]
                disp = "D3"
            elif funct.endswith("V"):
                base = funct[:-1]
                disp = "VV10"
            elif funct.endswith(")"):
                base = "-".join(funct.split("-")[0:-1])
                disp = "D3"
            else:
                base = funct
                disp = "None"

        avg = float(row[-1])

        if base in vac_mae:
            vac_mae[base][disp] = avg
        else:
            vac_mae[base] = {disp: avg}

with open(os.path.join(base_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]

        if funct in exclude:
            continue

        else:
            if funct.endswith("D"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D"):
                    base = base[:-1]
                disp = "D2"
            elif funct.endswith("D3"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D3"):
                    base = base[:-2]
                disp = "D3"
            elif funct.endswith("V"):
                base = funct[:-1]
                disp = "VV10"
            elif funct.endswith(")"):
                base = "-".join(funct.split("-")[0:-1])
                disp = "D3"
            else:
                base = funct
                disp = "None"

        avg = float(row[-1])

        if base in vac_rel:
            vac_rel[base][disp] = avg
        else:
            vac_rel[base] = {disp: avg}

with open(os.path.join(base_dir, "abserrs_IEF-PCM.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        if funct in exclude:
            continue

        else:
            if funct.endswith("D"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D"):
                    base = base[:-1]
                disp = "D2"
            elif funct.endswith("D3"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D3"):
                    base = base[:-2]
                disp = "D3"
            elif funct.endswith("V"):
                base = funct[:-1]
                disp = "VV10"
            elif funct.endswith(")"):
                base = "-".join(funct.split("-")[0:-1])
                disp = "D3"
            else:
                base = funct
                disp = "None"

        avg = float(row[-1])

        if base in pcm_mae:
            pcm_mae[base][disp] = avg
        else:
            pcm_mae[base] = {disp: avg}

with open(os.path.join(base_dir, "abserrs_rel_IEF-PCM.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        if funct in exclude:
            continue

        else:
            if funct.endswith("D"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D"):
                    base = base[:-1]
                disp = "D2"
            elif funct.endswith("D3"):
                if "-" in funct:
                    base = "-".join(funct.split("-")[0:-1])
                else:
                    base = funct

                if base.endswith("D3"):
                    base = base[:-2]
                disp = "D3"
            elif funct.endswith("V"):
                base = funct[:-1]
                disp = "VV10"
            elif funct.endswith(")"):
                base = "-".join(funct.split("-")[0:-1])
                disp = "D3"
            else:
                base = funct
                disp = "None"

        avg = float(row[-1])

        if base in pcm_rel:
            pcm_rel[base][disp] = avg
        else:
            pcm_rel[base] = {disp: avg}


fig, axs = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

for i, dset in enumerate([vac_mae, vac_rel, pcm_mae, pcm_rel]):
    ax = axs[i]

    if i % 2 == 0:
        ax.set_ylabel("MAE (eV)")
    else:
        ax.set_ylabel("MRAE (unitless)")

    x_labels = [x[0] for x in sorted(dset.items(), key=lambda x: x[0])]
    xs = np.arange(len(x_labels)) * 2

    for ii, disp in enumerate(["None", "D2", "D3", "VV10"]):
        this_ys = list()
        for x in x_labels:
            this_ys.append(dset[x].get(disp, 0.0))

        ax.bar(xs + (0.5 * ii - 0.75), this_ys, width=0.5, color=colors[disp], label=disp)

    if i == 3:
        ax.set_xticks(xs, x_labels, rotation=30, ha="right")

plt.tight_layout()
# plt.legend()

fig.savefig("sp_dispersion_effect.png", dpi=150)
# plt.show()