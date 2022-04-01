import csv
import os
import statistics
from monty.serialization import loadfn, dumpfn

import matplotlib.pyplot as plt
import numpy as np


def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


SMALL_SIZE = 12
MEDIUM_SIZE = 14
LARGE_SIZE = 18

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE,
       titlesize=LARGE_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels

base_dir = "../data"

methods = {
    "GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3",
            "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10"],
    "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS",
                 "TPSS-D3(BJ)", "MN12-L", "MN12-L-D3(BJ)", "B97M-rV"],
    "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP",
                   "CAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "wB97X",
                   "wB97XD", "wB97XD3", "wB97XV"],
    "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "M06-HF", "M08-SO", "M11",
                        "MN15", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)",
                        "SCAN0", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

g_mae_tzvpd = {x: dict() for x in methods}
h_mae_tzvpd = {x: dict() for x in methods}
s_mae_tzvpd = {x: dict() for x in methods}
e_mae_tzvpd = {x: dict() for x in methods}

for funct, avg in loadfn(os.path.join(base_dir,
                                      'g_mae_tzvpd_opts.json')).items():
    for group, functs in methods.items():
        if funct in functs:
            g_mae_tzvpd[group][funct] = avg
for funct, avg in loadfn(os.path.join(base_dir,
                                      'h_mae_tzvpd_opts.json')).items():
    for group, functs in methods.items():
        if funct in functs:
            h_mae_tzvpd[group][funct] = avg
for funct, avg in loadfn(os.path.join(base_dir,
                                      's_mae_tzvpd_opts.json')).items():
    for group, functs in methods.items():
        if funct in functs:
            s_mae_tzvpd[group][funct] = avg
for funct, avg in loadfn(os.path.join(base_dir,
                                      'elec_mae_tzvpd_opts.json')).items():
    for group, functs in methods.items():
        if funct in functs:
            e_mae_tzvpd[group][funct] = avg

fig, axs = plt.subplots(2, 2, figsize=(14, 6), sharex=True)

for i, dset in enumerate([e_mae_tzvpd, g_mae_tzvpd, h_mae_tzvpd, s_mae_tzvpd
                          ]):
    print(i)
    x = i // 2
    y = i % 2
    ax = axs[x][y]
    ax.set_ylabel("MAE (eV)")

    xs = ["GGA", "meta-GGA", "hybrid GGA", "hybrid meta-GGA"]
    avgs = list()
    lowlims = list()
    uplims = list()
    # for group in xs:
    #     avg = statistics.mean(dset[group].values())
    #     avgs.append(avg)
    #     group_sort = sorted(dset[group].items(), key=lambda x: x[1])
    #     print("\t min: {} ({}) max: {} ({}) avg: {}".format(group_sort[0][0],
    #                                                         group_sort[0][1],
    #                                                         group_sort[-1][0],
    #                                                         group_sort[-1][1],
    #                                                         avg))
    #     lowlims.append(abs(avg - group_sort[0][1]))
    #     uplims.append(abs(avg - group_sort[-1][1]))
    data = [sorted(np.array([v for k2, v in l.items()])) for k,l in dset.items()]
    parts = ax.violinplot(dataset=data, showmeans=True)
    colors = ["#ff595e", "#ffca3a", "#8ac926", "#1982c4"]
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    data2 = [sorted(np.random.normal(0, std, 100)) for std in range(1, 5)]
    print(data)
    print(data2)
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)

    # quartile1, medians, quartile3 = [], [], []
    # for d in data:
    #     q1, m, q3 = np.percentile(a=d, q=[25, 50, 75],
    #                                               axis=0)
    #     quartile1.append(q1)
    #     medians.append(m)
    #     quartile3.append(q3)
    whiskers = [
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data, quartile1, quartile3)]
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-',
              lw=1)

    set_axis_style(ax, labels=xs)

    #
    # ax.bar(x=xs, height=avgs, yerr=[lowlims, uplims],
    #        color=["#ff595e", "#ffca3a", "#8ac926", "#1982c4"])
    ax.set_xticklabels(xs, rotation=30, ha="right")

plt.tight_layout()
fig.savefig("average_tzvpd_thermochemistry_performance_opts.png", dpi=150)

plt.show()
