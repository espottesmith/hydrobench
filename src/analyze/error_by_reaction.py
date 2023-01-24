import os
import statistics

abs_file = "../data/results/new/abserrs_vacuum.csv"
rel_file = "../data/results/new/abserrs_rel_vacuum.csv"

methods = {"GGA": ["PBE", "PBE-D3(BJ)", "BLYP", "BLYP-D3(BJ)", "B97-D", "B97-D3", "mPW91", "mPW91-D3(BJ)", "VV10", "rVV10",],
           "meta-GGA": ["M06-L", "M06-L-D3(0)", "SCAN", "SCAN-D3(BJ)", "TPSS", "TPSS-D3(BJ)", "MN12-L", ("MN12-L-D3(BJ)"), "B97M-rV",],
           "hybrid GGA": ["PBE0", "PBE0-D3(BJ)", "LRC-wPBE", "LRC-wPBE-D3(BJ)", "LRC-wPBEh", "LRC-wPBEh-D3(BJ)", "B3LYP", "B3LYP-D3(BJ)", "CAM-B3LYP", "CAM-B3LYP-D3(0)", "rCAM-B3LYP", "rCAM-B3LYP-D3(0)", "mPW1PW91", "mPW1PW91-D3(BJ)", "HSE-HJS", "HSE-HJS-D3(BJ)", "wB97X", "wB97XD", "wB97XD3", "wB97XV",],
           "hybrid meta-GGA": ["M06-2X", "M06-2X-D3(0)", "wM06-D3", "M06-SX", "M06-SX-D3(BJ)", "M06-HF", "M06-HF-D3(0)", "M08-SO", "M08-SO-D3(0)", "M11", "M11-D3(0)", "revM11", "revM11-D3(0)", "MN15", "MN15-D3(0)", "BMK", "BMK-D3(BJ)", "TPSSh", "TPSSh-D3(BJ)", "SCAN0", "SCAN0-D3(BJ)", "mPWB1K", "mPWB1K-D3(BJ)", "wB97M-V"]}

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

        print(f"{reaction} & {statistics.mean(d['GGA']):.3f} & {statistics.mean(d['meta-GGA']):.3f} & {statistics.mean(d['hybrid GGA']):.3f} & {statistics.mean(d['hybrid meta-GGA']):.3f} & {statistics.mean(d['all']):.3f}\\\\\n\\hline")

print("\n\n\n")

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

        print(f"{reaction} & {statistics.mean(d['GGA']):.3f} & {statistics.mean(d['meta-GGA']):.3f} & {statistics.mean(d['hybrid GGA']):.3f} & {statistics.mean(d['hybrid meta-GGA']):.3f} & {statistics.mean(d['all']):.3f}\\\\\n\\hline")
