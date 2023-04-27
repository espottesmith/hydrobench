import os
import statistics

import numpy as np
import matplotlib.pyplot as plt


SMALL_SIZE = 16
MEDIUM_SIZE = 18
LARGE_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('title', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE, titlesize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels


abs_file = "../data/results/new/abserrs_vacuum.csv"
rel_file = "../data/results/new/abserrs_rel_vacuum.csv"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10",],
           "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", ("MN12-L-D3(BJ)"), "B97M-rV",],
           "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "LRC-wPBE", "LRC-wPBE-D3(BJ)", "LRC-wPBEh", "LRC-wPBEh-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "rCAM-B3LYP", "rCAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "HSE-HJS", "HSE-HJS-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV",],
           "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "wM06-D3", "M06-SX", "M06-SX-D3(BJ)", "M06-HF", "M06-HF-D3(0)", "M08-SO", "M08-SO-D3(0)", "M11", "M11-D3(0)", "revM11", "revM11-D3(0)", "MN15", "MN15-D3(0)", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "SCAN0-D3(BJ)", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

colors = {"GGA": "#D81B60", "meta-GGA": "#FFC107", "hybrid GGA": "#004D40", "hybrid meta-GGA": "#1E88E5", "all": "#2F4858"}

forward = list()
reverse = list()

fig, axs = plt.subplots(2, 2, figsize=(12, 10))

with open(abs_file, 'r') as file:
    file_contents = file.readlines()
    reactions = file_contents[0].split(",")[1:-1]

    data = dict()

    for reaction in reactions:
        data[reaction] = {"all": list()}
        for funct_class in methods:
            data[reaction][funct_class] = list()

    for line in file_contents[1:]:
        line_contents = line.split(",")
        functional = line_contents[0]

        if functional in ["r2SCAN3c"]:
            continue

        this_class = None
        for funct_class, functs in methods.items():
            if functional in functs:
                this_class = funct_class
                break
        if this_class is None:
            print("CANNOT FIND FUNCTIONAL CLASS", functional)
            continue

        errors = list()
        for x in line_contents[1:-1]:
            try:
                errors.append(float(x))
            except:
                errors.append(None)

        for ii, reaction in enumerate(reactions):
            error = errors[ii]
            data[reaction]["all"].append(error)
            data[reaction][this_class].append(error)

    for reaction, d in data.items():
        if None in d["all"]:
            continue

        if any([x in reaction for x in ["amide11", "amide12", "enamine", "amide22", "furan3"]]):
            continue

        if "forward" in reaction:
            forward.append((reaction, {
                "GGA": statistics.mean(d['GGA']),
                "meta-GGA": statistics.mean(d['meta-GGA']),
                "hybrid GGA": statistics.mean(d['hybrid GGA']),
                "hybrid meta-GGA": statistics.mean(d["hybrid meta-GGA"]),
                "all": statistics.mean(d['all'])
            }))
        else:
            reverse.append((reaction, {
                "GGA": statistics.mean(d['GGA']),
                "meta-GGA": statistics.mean(d['meta-GGA']),
                "hybrid GGA": statistics.mean(d['hybrid GGA']),
                "hybrid meta-GGA": statistics.mean(d["hybrid meta-GGA"]),
                "all": statistics.mean(d['all'])
            }))

        # print(f"{reaction} & {statistics.mean(d['GGA']):.3f} & {statistics.mean(d['meta-GGA']):.3f} & {statistics.mean(d['hybrid GGA']):.3f} & {statistics.mean(d['hybrid meta-GGA']):.3f} & {statistics.mean(d['all']):.3f}\\\\\n\\hline")

# forward first
print("FORWARD MAE")
dsets = ["GGA", "meta-GGA", "hybrid GGA", "hybrid meta-GGA", "all"]

# Print averages
for dset in dsets:
    print(f"{dset} FORWARD: {statistics.mean([x[1][dset] for x in forward])} REVERSE: {statistics.mean([x[1][dset] for x in reverse])}")

box_values = list()
for dset in dsets:
    print("\t", dset)
    d = np.array(sorted([r[1][dset] for r in forward]))
    box_values.append(d)
    q1 = np.quantile(d, 0.25)
    q3 = np.quantile(d, 0.75)
    med = np.median(d)
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)
    lower_bound = q1 - (1.5 * iqr)

    for n, v in forward:
        if v[dset] <= lower_bound or v[dset] >= upper_bound:
            print("\t\t", n, v)

bp = axs[0][0].boxplot(box_values,
                    #    labels=dsets,
                       patch_artist=True,
                       widths=0.6)
axs[0][0].set_ylabel("MAE (eV)")
axs[0][0].set_xticks([])
axs[0][0].set_ylim(0, 0.45)

for patch, color in zip(bp['boxes'], [colors[x] for x in dsets]):
    patch.set_facecolor(color)

for median in bp['medians']:
    median.set(color='black')

box_values = list()
for dset in dsets:
    print("\t", dset)
    d = np.array(sorted([r[1][dset] for r in reverse]))
    box_values.append(d)
    q1 = np.quantile(d, 0.25)
    q3 = np.quantile(d, 0.75)
    med = np.median(d)
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)
    lower_bound = q1 - (1.5 * iqr)

    for n, v in reverse:
        if v[dset] <= lower_bound or v[dset] >= upper_bound:
            print("\t\t", n, v)

bp = axs[0][1].boxplot(box_values,
                    #    labels=dsets,
                       patch_artist=True,
                       widths=0.6)
axs[0][1].set_ylabel("MAE (eV)")
axs[0][1].set_xticks([])
axs[0][1].set_ylim(0, 0.45)

for patch, color in zip(bp['boxes'], [colors[x] for x in dsets]):
    patch.set_facecolor(color)

for median in bp['medians']:
    median.set(color='black')

print("\n\n\n")

forward = list()
reverse = list()


with open(rel_file, 'r') as file:
    file_contents = file.readlines()
    reactions = file_contents[0].split(",")[1:-1]

    data = dict()

    for reaction in reactions:
        data[reaction] = {"all": list()}
        for funct_class in methods:
            data[reaction][funct_class] = list()

    for line in file_contents[1:]:
        line_contents = line.split(",")
        functional = line_contents[0]

        if functional in ["r2SCAN3c"]:
            continue

        this_class = None
        for funct_class, functs in methods.items():
            if functional in functs:
                this_class = funct_class
                break
        if this_class is None:
            print("CANNOT FIND FUNCTIONAL CLASS", functional)
            continue

        errors = list()
        for x in line_contents[1:-1]:
            try:
                errors.append(float(x))
            except:
                errors.append(None)

        for ii, reaction in enumerate(reactions):
            error = errors[ii]
            data[reaction]["all"].append(error)
            data[reaction][this_class].append(error)

    for reaction, d in data.items():
        if None in d["all"]:
            continue

        if "forward" in reaction:
            forward.append((reaction, {
                "GGA": statistics.mean(d['GGA']),
                "meta-GGA": statistics.mean(d['meta-GGA']),
                "hybrid GGA": statistics.mean(d['hybrid GGA']),
                "hybrid meta-GGA": statistics.mean(d["hybrid meta-GGA"]),
                "all": statistics.mean(d['all'])
            }))
        else:
            reverse.append((reaction, {
                "GGA": statistics.mean(d['GGA']),
                "meta-GGA": statistics.mean(d['meta-GGA']),
                "hybrid GGA": statistics.mean(d['hybrid GGA']),
                "hybrid meta-GGA": statistics.mean(d["hybrid meta-GGA"]),
                "all": statistics.mean(d['all'])
            }))

        # print(f"{reaction} & {statistics.mean(d['GGA']):.3f} & {statistics.mean(d['meta-GGA']):.3f} & {statistics.mean(d['hybrid GGA']):.3f} & {statistics.mean(d['hybrid meta-GGA']):.3f} & {statistics.mean(d['all']):.3f}\\\\\n\\hline")

# forward first
print("FORWARD MAE")
dsets = ["GGA", "meta-GGA", "hybrid GGA", "hybrid meta-GGA", "all"]
box_values = list()
for dset in dsets:
    print("\t", dset)
    d = np.array(sorted([r[1][dset] for r in forward]))
    box_values.append(d)
    q1 = np.quantile(d, 0.25)
    q3 = np.quantile(d, 0.75)
    med = np.median(d)
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)
    lower_bound = q1 - (1.5 * iqr)

    for n, v in forward:
        if v[dset] <= lower_bound or v[dset] >= upper_bound:
            print("\t\t", n, v)

bp = axs[1][0].boxplot(box_values,
                    #    labels=dsets,
                       patch_artist=True,
                       widths=0.6)
axs[1][0].set_ylabel("MRAE (unitless)")
axs[1][0].set_xticks([])

for patch, color in zip(bp['boxes'], [colors[x] for x in dsets]):
    patch.set_facecolor(color)

for median in bp['medians']:
    median.set(color='black')

box_values = list()
for dset in dsets:
    print("\t", dset)
    d = np.array(sorted([r[1][dset] for r in reverse]))
    box_values.append(d)
    q1 = np.quantile(d, 0.25)
    q3 = np.quantile(d, 0.75)
    med = np.median(d)
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)
    lower_bound = q1 - (1.5 * iqr)

    for n, v in reverse:
        if v[dset] <= lower_bound or v[dset] >= upper_bound:
            print("\t\t", n, v)

bp = axs[1][1].boxplot(box_values,
                    #    labels=dsets,
                       patch_artist=True,
                       widths=0.6)
axs[1][1].set_ylabel("MRAE (unitless)")
axs[1][1].set_xticks([])

for patch, color in zip(bp['boxes'], [colors[x] for x in dsets]):
    patch.set_facecolor(color)

for median in bp['medians']:
    median.set(color='black')

plt.tight_layout(pad=3)
fig.savefig("forward_vs_reverse_v3.png", dpi=300)

plt.show()