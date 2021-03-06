import os
import statistics
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sbn

from monty.serialization import loadfn

SMALL_SIZE = 12
MEDIUM_SIZE = 14
LARGE_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

dft_methods = ["SCAN0", "mPWB1K", ("mPWB1K", "D3(BJ)"), "wB97M-V"]

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10"],
            "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", "MN12-L-D3(BJ)", "B97M-rV"],
            "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", r"$\omega$B97X", r"$\omega$B97X-D", r"$\omega$B97X-D3", r"$\omega$B97X-V"],
            "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "M06-HF", "M08-SO", "M11", "MN15", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "mPWB1K", "mPWB1K-D3(BJ)", r"$\omega$B97M-V"]}

colors = {"GGA": "#ff595e",
          "meta-GGA": "#ffca3a",
          "hybrid GGA": "#8ac926",
          "hybrid meta-GGA": "#1982c4"}

replacements = {"pbe": "PBE",
                "blyp": "BLYP",
                "b97-d": "B97-D",
                "b97-d3": "B97-D3",
                "mpw91": "mPW91",
                "vv10": "VV10",
                "rvv10": "rVV10",
                "m06-l": "M06-L",
                "scan": "SCAN",
                "tpss": "TPSS",
                "mn12-l": "MN12-L",
                "b97m-rv": "B97M-rV",
                "pbe0": "PBE0",
                "b3lyp": "B3LYP",
                "cam-b3lyp": "CAM-B3LYP",
                "mpw1pw91": "mPW1PW91",
                "wb97x": r"$\omega$B97X",
                "wb97xd": r"$\omega$B97X-D",
                "wb97xd3": r"$\omega$B97X-D3",
                "wb97xv": r"$\omega$B97X-V",
                "m06-2x": "M06-2X",
                "m06-hf": "M06-HF",
                "m08-so": "M08-SO",
                "m11": "M11",
                "mn15": "MN15",
                "bmk": "BMK",
                "tpssh": "TPSSh",
                "scan0": "SCAN0",
                "mpwb1k": "mPWB1K",
                "wb97m-v": r"$\omega$B97M-V"
                }

sp_data = loadfn("/Users/ewcss/data/ssbt/20220204_sp_benchmark/dft_sp.json")

walltimes = defaultdict(list)
cputimes = defaultdict(list)

for sp in sp_data:
    method = sp["orig"]["rem"]["method"]
    method = replacements.get(method, method)
    if "dft_d" in sp["orig"]["rem"]:
        dft_d = sp["orig"]["rem"]["dft_d"]
        if dft_d.upper() == "D3_BJ":
            method += "-D3(BJ)"
        elif dft_d.upper() == "D3_ZERO":
            method += "-D3(0)"
    walltime = sp["walltime"]
    cputime = sp["walltime"]
    walltimes[method].append(walltime)
    cputimes[method].append(cputime)

wall = dict()
cpu = dict()

for m, w in walltimes.items():
    wall[m] = statistics.mean(w)
for m, c in cputimes.items():
    cpu[m] = statistics.mean(c)


min_wall = min(wall.values())
min_cpu = min(cpu.values())

for m, w in wall.items():
    wall[m] = wall[m] / min_wall
for m, c in cpu.items():
    cpu[m] = cpu[m] / min_cpu

fig, axs = plt.subplots(figsize=(11, 12))

sorted_wall = sorted(wall.items(), key=lambda x: x[1])
sorted_cpu = sorted(cpu.items(), key=lambda x: x[1])

c = list()
for h in sorted_wall:
    match = False
    for g, functs in methods.items():
        if h[0] in functs:
            match = True
            c.append(colors[g])
            break
    if not match:
        print("problem", h[0])

xs =  [h[0] for h in sorted_wall]
ypos = np.arange(len(xs))
axs.barh(np.arange(len(xs)), [h[1] for h in sorted_wall], color=c, align="center")
# axs.barh(y_pos, performance, xerr=error, align='center')
axs.set_yticks(ypos, labels=xs)
axs.invert_yaxis()
plt.title("Cost of Single-Point Energy Calculations")
plt.xlabel("Relative Walltime")
# plt.ylabel("Functional")
# sbn.barplot(x=[h[1] for h in sorted_cpu], y=[h[0] for h in sorted_cpu], ax=axs[1])
plt.tight_layout()

plt.show()
fig.savefig("relative_cost.png", dpi=150)