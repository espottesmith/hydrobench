import os
import re
import copy

import numpy as np

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

hf = re.compile(r"\s*SCF\s+energy in the final basis set\s+=\s+([\-\.0-9]+)")
mp2_corr = re.compile(r"\s*Total\s+MP2\s+correlation\s+energy\s+=\s+([\-\.0-9]+)")
mp2 = re.compile(r"\s+MP2\s+total\s+energy\s+=\s+([\-\.0-9]+)")

alpha_svp_tzvp = 10.390
alpha_tzvpp_qzvpp = 7.880

beta_svp_tzvp = 2.400
beta_tzvpp_qzvpp = 2.970

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)

all_data = dict()

for base_dir in ["/projects/silimorphous/ewcss/ssbt/20220209_mp2/block_2022-02-09-22-24-40-918618"]:
    for direct in os.listdir(base_dir):
        files = os.listdir(os.path.join(base_dir, direct))

        # Nothing ran
        if not any(["FW.json" in x for x in files]):
            continue
        # Error
        if "error.1.tar.gz" in files:
            continue

        if "FW.json" in files:
            meta = loadfn(os.path.join(base_dir, direct, "FW.json"))
        else:
            meta = loadfn(os.path.join(base_dir, direct, "FW.json.gz"))

        name = meta["name"]

        if "SMD" in name:
            continue

        all_data[name] = dict()

        if "mol.qout" in files:
            filename = "mol.qout"
        else:
            filename = "mol.qout.gz"

        with zopen(os.path.join(base_dir, direct, filename), mode="rt", encoding="ISO-8859-1") as f:
            text = f.read()

            thishf = hf.search(text)
            if thishf is not None:
                all_data[name]["hf"] = float(thishf.group(1))

            thismp2_corr = mp2_corr.search(text)
            if thismp2_corr is not None:
                all_data[name]["mp2_corr"] = float(thismp2_corr.group(1))

            thismp2 = mp2.search(text)
            if thismp2 is not None:
                all_data[name]["mp2"] = float(thismp2.group(1))

dumpfn(all_data, "all_mp2_data.json")

groups = dict()

for n, data in all_data.items():

    name = copy.deepcopy(n)

    for a, b in replacements.items():
        name = name.replace(a, b)

    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    try:
        basis = contents[4]
    except IndexError:
        print(n)
        continue

    if groupname not in groups:
        groups[groupname] = {basis: data}
    else:
        groups[groupname][basis] = data

dumpfn(groups, "grouped_mp2_data.json")

full_corrected_data = dict()

for groupname, groupdata in groups.items():
    try:
        hf = correct_hf(groupdata["def2-TZVPP"]["hf"],
                        groupdata["def2-QZVPP"]["hf"],
                        3,
                        4,
                        alpha=alpha_tzvpp_qzvpp)

        mp2 = correct_wft(groupdata["def2-TZVPP"]["mp2_corr"],
                           groupdata["def2-QZVPP"]["mp2_corr"],
                           3,
                           4,
                           beta=beta_tzvpp_qzvpp)

        total_e = hf + mp2
        full_corrected_data[groupname] = total_e
    except KeyError:
        print("PROBLEM CORRECTED FULL", groupname)
        continue

dumpfn(full_corrected_data, "mp2_corrected_full.json")