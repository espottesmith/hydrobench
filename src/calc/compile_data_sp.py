import copy
import os
import statistics

import numpy as np

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from atomate.qchem.database import QChemCalcDb


def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)


db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

base_dir = "/Users/ewcss/data/ssbt/20220204_sp_benchmark"

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

dft_data = [e for e in db.db["tasks"].find({"tags.set": "20211016_sp_dft"}) if "SMD" not in e["task_label"]]
cc_data = [e for e in db.db["tasks"].find({"tags.set": "20211016_sp_cc"}) if "SMD" not in e["task_label"]]

for calc in dft_data + cc_data:
    label = copy.deepcopy(calc["task_label"])
    for a, b in replacements.items():
        label = label.replace(a, b)

    calc["task_label"] = label

all_cc_data = {
    e["task_label"]: {
        "hf": e["calcs_reversed"][0]["hf_scf_energy"],
        "ccsd_corr": e["calcs_reversed"][0]["ccsd_correlation_energy"],
        "ccsd_total": e["calcs_reversed"][0]["ccsd_total_energy"],
        "ccsdt_corr": e["calcs_reversed"][0].get("ccsd(t)_total_energy")
    } for e in cc_data if "SMD" not in e["task_label"]}

dumpfn(dft_data, os.path.join(base_dir, "dft_sp.json"))
dumpfn(cc_data, os.path.join(base_dir, "cc_sp.json"))

alpha_svp_tzvp = 10.390
alpha_tzvpp_qzvpp = 7.880

beta_svp_tzvp = 2.400
beta_tzvpp_qzvpp = 2.970

dumpfn(all_cc_data, os.path.join(base_dir, "cc_data_parsed.json"))

groups = dict()

for name, data in all_cc_data.items():
    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    meth_basis = "_".join([contents[3], contents[4]])

    if groupname not in groups:
        groups[groupname] = {meth_basis: data}
    else:
        groups[groupname][meth_basis] = data

dumpfn(groups, os.path.join(base_dir, "cc_data_grouped.json"))

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
        print("PROBLEM FULLY CORRECTED", groupname)
        continue

dumpfn(full_corrected_data, os.path.join(base_dir, "cc_full_corrected.json"))

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

dumpfn(corrected_small_data, os.path.join(base_dir, "cc_corrected_small.json"))

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

dumpfn(corrected_data, os.path.join(base_dir, "cc_corrected_t.json"))