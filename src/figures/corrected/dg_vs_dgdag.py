import csv
import os
import difflib
import statistics

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

SMALL_SIZE = 16
MEDIUM_SIZE = 18
LARGE_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

bar_dir = "../../data/results/corrected/"
dg_dir = "../../data/results/dg/"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10",],
           "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", ("MN12-L-D3(BJ)"), "B97M-rV",],
           "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "LRC-wPBE", "LRC-wPBE-D3(BJ)", "LRC-wPBEh", "LRC-wPBEh-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "rCAM-B3LYP", "rCAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "HSE-HJS", "HSE-HJS-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV",],
           "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "wM06-D3", "M06-SX", "M06-SX-D3(BJ)", "M06-HF", "M06-HF-D3(0)", "M08-SO", "M08-SO-D3(0)", "M11", "M11-D3(0)", "revM11", "revM11-D3(0)", "MN15", "MN15-D3(0)", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "SCAN0-D3(BJ)", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

colors = {"GGA": "#D81B60", "meta-GGA": "#FFC107", "hybrid GGA": "#004D40", "hybrid meta-GGA": "#1E88E5"}

functs = dict()
functs_rel = dict()

with open(os.path.join(bar_dir, "abserrs_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        functs[funct] = {"barrier": avg}

with open(os.path.join(dg_dir, "abserrs_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        functs[funct]["dg"] = avg

with open(os.path.join(bar_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        functs_rel[funct] = {"barrier": avg}

with open(os.path.join(dg_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])

        functs_rel[funct]["dg"] = avg

fig, ax = plt.subplots(1, 2, figsize=(12, 6))

xs0 = np.linspace(0.00, 0.25, 100)
ax[0].plot(xs0, xs0, '--k')
ax[0].set_xlim(0, 0.25)
ax[0].set_ylim(0, 0.25)
ax[0].fill_between(xs0, xs0 * 1.5, xs0 * 3, color="#DCDCDC")
ax[0].fill_between(xs0, xs0 / 3, xs0 / 1.5, color="#DCDCDC")
ax[0].fill_between(xs0, xs0 * 3, 1.0, color="#C0C0C0")
ax[0].fill_between(xs0, -0.1, xs0 / 3, color="#C0C0C0")


for group, functionals in methods.items():
    print(group)
    color = colors.get(group, "#000000")
    xs = list()
    ys = list()
    for f in functionals:
        data = functs.get(f)
        print("\t", f, data)
        xs.append(data["dg"])
        ys.append(data["barrier"])

    ax[0].scatter(xs, ys, color=color, alpha=0.8, edgecolors="black", s=80, label=group)

ax[0].set_xlabel(r"MAE$_{\Delta E}$ (eV)")
ax[0].set_ylabel(r"MAE$_{\Delta E^‡}$ (eV)")
ax[0].set_aspect("equal", "box")

print()
print()

xs1 = np.linspace(0.00, 0.9, 100)
ax[1].plot(xs1, xs1, '--k')
ax[1].set_xlim(0, 0.9)
ax[1].set_ylim(0, 0.9)
ax[1].fill_between(xs1, xs1 * 1.5, xs1 * 3, color="#DCDCDC")
ax[1].fill_between(xs1, xs1 / 3, xs1 / 1.5, color="#DCDCDC")
ax[1].fill_between(xs1, xs1 * 3, 1.0, color="#C0C0C0")
ax[1].fill_between(xs1, -0.1, xs1 / 3, color="#C0C0C0")


for group, functionals in methods.items():
    print(group)
    color = colors.get(group, "#000000")
    xs = list()
    ys = list()
    for f in functionals:
        data = functs_rel.get(f)
        print("\t", f, data)
        xs.append(data["dg"])
        ys.append(data["barrier"])

    ax[1].scatter(xs, ys, color=color, alpha=0.8, edgecolors="black", s=80, label=group)

ax[1].set_xlabel(r"MRAE$_{\Delta E}$ (unitless)")
ax[1].set_ylabel(r"MRAE$_{\Delta E^‡}$ (unitless)")
ax[1].set_aspect("equal", "box")

ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# ax[1].legend()
plt.tight_layout()
fig.savefig("dg_vs_dgdag_eq.png", dpi=300)
plt.show()
