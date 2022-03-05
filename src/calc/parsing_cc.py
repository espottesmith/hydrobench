import os
import re

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

def correct_hf(e1, e2, n1, n2, alpha=7.880):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta=2.970):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)


all_data = loadfn("/Users/ewcss/data/ssbt/20211023_cc/all_cc_data.json")

groups = dict()

for name, data in all_data.items():
    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    meth_basis = "_".join([contents[3], contents[4]])

    if groupname not in groups:
        groups[groupname] = {meth_basis: data}
    else:
        groups[groupname][meth_basis] = data

dumpfn(groups, "/Users/ewcss/data/ssbt/20220201_cc/grouped_cc_data.json")

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
        continue

dumpfn(full_corrected_data, "/Users/ewcss/data/ssbt/20220201_cc/full_corrected_cc_data.json")

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
        continue

dumpfn(corrected_data, "/Users/ewcss/data/ssbt/20220201_cc/corrected_t_cc_data.json")

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
        continue

dumpfn(corrected_small_data, "/Users/ewcss/data/ssbt/20220201_cc/corrected_small_cc_data.json")

uncorrected_data = dict()
for groupname, groupdata in groups.items():
    try:
        total_ccsd_e = groupdata["CCSD_def2-TZVPP"]["ccsd_total"]
        uncorrected_data[groupname] = total_ccsd_e
    except KeyError:
        continue

dumpfn(uncorrected_data, "/Users/ewcss/data/ssbt/20220201_cc/uncorrected_cc_data.json")