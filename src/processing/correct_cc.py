import re
import numpy as np

from monty.serialization import loadfn, dumpfn

alpha_svp_tzvp = 10.390
alpha_tzvpp_qzvpp = 7.880

beta_svp_tzvp = 2.400
beta_tzvpp_qzvpp = 2.970

def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)


all_data = loadfn("../data/reference/all_cc_data.json")
corrected_data = loadfn("../data/reference/corrected/20230622_dumped_cc.json")

for k, v in corrected_data.items():
    new_k = k.replace("realts", "ts")
    if new_k in all_data:
        print("REPLACING", new_k)
        all_data[new_k] = v

groups = dict()

for name, data in all_data.items():
    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    meth_basis = "_".join([contents[3], contents[4]])

    if groupname not in groups:
        groups[groupname] = {meth_basis: data}
    else:
        groups[groupname][meth_basis] = data

dumpfn(groups, "../data/reference/corrected/grouped_cc_data.json")

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

dumpfn(full_corrected_data, "../data/reference/corrected/full_corrected_cc_data.json")

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

dumpfn(corrected_data, "../data/reference/corrected/corrected_t_cc_data.json")

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

dumpfn(corrected_small_data, "../data/reference/corrected/corrected_small_cc_data.json")

uncorrected_data = dict()
for groupname, groupdata in groups.items():
    try:
        total_ccsd_e = groupdata["CCSD_def2-TZVPP"]["ccsd_total"]
        uncorrected_data[groupname] = total_ccsd_e
    except KeyError:
        continue

dumpfn(uncorrected_data, "../data/reference/corrected/uncorrected_cc_data.json")