import os
import re
import copy

import numpy as np

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

hf = re.compile(r"\s+SCF energy\s+=\s+([\-\.0-9]+)")
mp2 = re.compile(r"\s+MP2 energy\s+=\s+([\-\.0-9]+)")
ccsd_corr = re.compile(r"\s+CCSD correlation energy\s+=\s+([\-\.0-9]+)")
ccsd_total = re.compile(r"\s+CCSD total energy\s+=\s+([\-\.0-9]+)")
ccsdt_corr = re.compile(r"\s+CCSD\(T\) correlation energy\s+=\s+([\-\.0-9]+)")
ccsdt_total = re.compile(r"\s+CCSD\(T\) total energy\s+=\s+([\-\.0-9]+)")

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

for base_dir in ["/projects/silimorphous/ewcss/ssbt/20211016_sp_cc/block_2021-10-18-14-42-59-164123",
                "/projects/silimorphous/ewcss/ssbt/20211016_sp_cc/block_2022-01-27-22-38-26-594047",
                 "/projects/silimorphous/ewcss/ssbt/20211021_sp_cc/block_2021-10-21-20-34-57-567627",
                 ]:
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

            thismp2 = mp2.search(text)
            if thismp2 is not None:
                all_data[name]["mp2"] = float(thismp2.group(1))

            thisccsd_corr = ccsd_corr.search(text)
            if thisccsd_corr is not None:
                all_data[name]["ccsd_corr"] = float(thisccsd_corr.group(1))

            thisccsd_total = ccsd_total.search(text)
            if thisccsd_total is not None:
                all_data[name]["ccsd_total"] = float(thisccsd_total.group(1))

            thisccsdt_corr = ccsdt_corr.search(text)
            if thisccsdt_corr is not None:
                all_data[name]["ccsdt_corr"] = float(thisccsdt_corr.group(1))

            thisccsdt_total = ccsdt_total.search(text)
            if thisccsdt_total is not None:
                all_data[name]["ccsdt_total"] = float(thisccsdt_total.group(1))

dumpfn(all_data, "all_cc_data.json")

groups = dict()

for n, data in all_data.items():

    name = copy.deepcopy(n)

    for a, b in replacements.items():
        name = name.replace(a, b)

    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    meth_basis = "_".join([contents[3], contents[4]])

    if groupname not in groups:
        groups[groupname] = {meth_basis: data}
    else:
        groups[groupname][meth_basis] = data

dumpfn(groups, "grouped_cc_data.json")

full_corrected_data = dict()

for groupname, groupdata in groups.items():
    try:
        hf = correct_hf(groupdata["CCSD_def2-TZVPP"]["hf"],
                        groupdata["CCSD_def2-QZVPP"]["hf"],
                        3,
                        4,
                        alpha=alpha_tzvpp_qzvpp)

        ccsd = correct_wft(groupdata["CCSD_def2-TZVPP"]["ccsd_corr"],
                           groupdata["CCSD_def2-QZVPP"]["ccsd_corr"],
                           3,
                           4,
                           beta=beta_tzvpp_qzvpp)

        paren_t = correct_wft(groupdata["CCSD(T)_def2-SVP"]["ccsdt_corr"],
                              groupdata["CCSD(T)_def2-TZVP"]["ccsdt_corr"],
                              2,
                              3,
                              beta=beta_svp_tzvp)
        total_e = hf + ccsd + paren_t
        full_corrected_data[groupname] = total_e
    except KeyError:
        print("PROBLEM CORRECTED SMALL", groupname)
        continue

dumpfn(full_corrected_data, "cc_corrected_full.json")

# Not currently correcting HF or CCSD
# Only applying an extrapolation to (T), due to issues with QZVPP calculations
corrected_data = dict()
for groupname, groupdata in groups.items():
    try:
        hf_ccsd = groupdata["CCSD_def2-TZVPP"]["ccsd_total"]
        paren_t = correct_wft(groupdata["CCSD(T)_def2-SVP"]["ccsdt_corr"],
                              groupdata["CCSD(T)_def2-TZVP"]["ccsdt_corr"],
                              2,
                              3,
                              beta=beta_svp_tzvp)
        total_e = hf_ccsd + paren_t
        corrected_data[groupname] = total_e
    except KeyError:
        print("PROBLEM CORRECTED T", groupname)
        continue

dumpfn(corrected_data, "cc_corrected_t.json")

uncorrected_data = dict()
for groupname, groupdata in groups.items():
    try:
        total_ccsd_e = groupdata["CCSD_def2-TZVPP"]["ccsd_total"]
        uncorrected_data[groupname] = total_ccsd_e
    except KeyError:
        print("PROBLEM UNCORRECTED", groupname)
        continue

dumpfn(uncorrected_data, "ccsd_uncorrected.json")

corrected_small_data = dict()
for groupname, groupdata in groups.items():
    try:
        hf = correct_hf(groupdata["CCSD(T)_def2-SVP"]["hf"],
                        groupdata["CCSD(T)_def2-TZVP"]["hf"],
                        3,
                        4,
                        alpha=alpha_svp_tzvp)

        ccsd = correct_wft(groupdata["CCSD(T)_def2-SVP"]["ccsd_corr"],
                           groupdata["CCSD(T)_def2-TZVP"]["ccsd_corr"],
                           3,
                           4,
                           beta=beta_svp_tzvp)

        paren_t = correct_wft(groupdata["CCSD(T)_def2-SVP"]["ccsdt_corr"],
                              groupdata["CCSD(T)_def2-TZVP"]["ccsdt_corr"],
                              2,
                              3,
                              beta=beta_svp_tzvp)
        total_e = hf + ccsd + paren_t
        corrected_small_data[groupname] = total_e
    except KeyError:
        print("PROBLEM CORRECTED SMALL", groupname)
        continue

dumpfn(corrected_small_data, "cc_corrected_small.json")