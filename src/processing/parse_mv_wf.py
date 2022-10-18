import os
import re

import numpy as np

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.database import QChemCalcDb

from .utils import (correct_hf, correct_wft,
                    alpha_tzvpp_qzvpp,
                    beta_tzvpp_qzvpp)

db = QChemCalcDb.from_db_file("db.json")

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

# For MP2 calculations
hf_mp2 = re.compile(r"\s*SCF\s+energy in the final basis set\s+=\s+([\-\.0-9]+)")
mp2_mp2 = re.compile(r"\s*Total\s+MP2\s+correlation\s+energy\s+=\s+([\-\.0-9]+)")
mp2_total = re.compile(r"\s+MP2\s+total\s+energy\s+=\s+([\-\.0-9]+)")

hf = re.compile(r"\s+SCF energy\s+=\s+([\-\.0-9]+)")
mp2 = re.compile(r"\s+MP2(?:\s+total)?\s+energy\s+=\s+([\-\.0-9]+)")
ccsd_corr = re.compile(r"\s+CCSD correlation energy\s+=\s+([\-\.0-9]+)")
ccsd_total = re.compile(r"\s+CCSD total energy\s+=\s+([\-\.0-9]+)")
ccsdt_corr = re.compile(r"\s+CCSD\(T\) correlation energy\s+=\s+([\-\.0-9]+)")
ccsdt_total = re.compile(r"\s+CCSD\(T\) total energy\s+=\s+([\-\.0-9]+)")


base_dir = "20220309_mv_sp_wf_dir"

grouped_data = dict()

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

    for a, b in replacements.items():
        name = name.replace(a, b)

    method = meta["spec"]["_tasks"][0]["qchem_input_params"]["overwrite_inputs"]["rem"]["method"]
    basis = meta["spec"]["_tasks"][0]["qchem_input_params"]["basis_set"]
    if method.lower() == "mp2":
        is_mp2 = True
    else:
        is_mp2 = False

    if "SMD" in name:
        continue

    if "mol.qout" in files:
        filename = "mol.qout"
    else:
        filename = "mol.qout.gz"

    qcout = QCOutput(os.path.join(base_dir, direct, filename))
    mol = qcout.data["initial_molecule"]
    calc_mg = MoleculeGraph.with_local_env_strategy(
        mol,
        OpenBabelNN()
    )

    this_rxn = name.split("_")[0]
    this_state = name.split("_")[1]
    this_solv = name.split("_")[-1]

    if this_rxn not in grouped_data:
        grouped_data[this_rxn] = {solv : {"ts": dict(),
                                  "rct": dict(),
                                  "pro": dict()} for solv in ["vacuum", "IEF-PCM"]}

    with zopen(os.path.join(base_dir, direct, filename), mode="rt", encoding="ISO-8859-1") as f:
        text = f.read()

        if is_mp2:
            thishf = hf_mp2.search(text)
            if thishf is not None:
                grouped_data[this_rxn][this_solv][this_state]["hf_" + basis] = float(thishf.group(1))

            thismp2_corr = mp2_mp2.search(text)
            if thismp2_corr is not None:
                grouped_data[this_rxn][this_solv][this_state]["mp2_corr_" + basis] = float(thismp2_corr.group(1))

            thismp2 = mp2_total.search(text)
            if thismp2 is not None:
                grouped_data[this_rxn][this_solv][this_state]["mp2_total_" + basis] = float(thismp2.group(1))

        else:
            thishf = hf.search(text)
            if thishf is not None:
                grouped_data[this_rxn][this_solv][this_state]["hf_" + basis] = float(thishf.group(1))

            thismp2 = mp2.search(text)
            if thismp2 is not None:
                grouped_data[this_rxn][this_solv][this_state]["ccsdt_mp2_" + basis] = float(thismp2.group(1))

            thisccsd_corr = ccsd_corr.search(text)
            if thisccsd_corr is not None:
                grouped_data[this_rxn][this_solv][this_state]["ccsd_corr_" + basis] = float(thisccsd_corr.group(1))

            thisccsd_total = ccsd_total.search(text)
            if thisccsd_total is not None:
                grouped_data[this_rxn][this_solv][this_state]["ccsd_total_" + basis] = float(thisccsd_total.group(1))

            thisccsdt_corr = ccsdt_corr.search(text)
            if thisccsdt_corr is not None:
                grouped_data[this_rxn][this_solv][this_state]["ccsdt_corr_" + basis] = float(thisccsdt_corr.group(1))

            thisccsdt_total = ccsdt_total.search(text)
            if thisccsdt_total is not None:
                grouped_data[this_rxn][this_solv][this_state]["ccsdt_total_" + basis] = float(thisccsdt_total.group(1))

full_corrected_data = dict()

for groupname, groupdata in grouped_data.items():
    for solvname, solvdata in groupdata.items():
        for statename, statedata in solvdata.items():
            # print(groupname, statename, statedata.keys())
            hf = correct_hf(statedata["hf_def2-TZVPP"],
                            statedata["hf_def2-QZVPP"],
                            3,
                            4,
                            alpha=alpha_tzvpp_qzvpp)

            mp2 = correct_wft(statedata["mp2_corr_def2-TZVPP"],
                              statedata["mp2_corr_def2-QZVPP"],
                              3,
                              4,
                              beta=beta_tzvpp_qzvpp)

            total_e = hf + mp2

            try:
                diff = statedata["ccsdt_total_def2-TZVP"] - statedata["ccsdt_mp2_def2-TZVP"]
            except:
                # One case where data wasn't parsed properly
                # Manually extracted relevant energies from output files
                print("PROBLEM", groupname, solvname, statename)
                diff = -416.61050695 - -416.50353665

            grouped_data[groupname][solvname][statename]["total_corrected"] = total_e + diff

dumpfn(grouped_data, "../data/reference/mvopt_svpd_reference_data.json")