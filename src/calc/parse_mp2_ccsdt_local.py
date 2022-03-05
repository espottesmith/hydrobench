import os
import re
import copy

import numpy as np

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.database import QChemCalcDb

base_dir_mols = "/Users/ewcss/data/ssbt/for_sp"

baseline_mg = {
    "amide22": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "amide_2_2_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "amide_2_2_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "amide_2_2_pro_-1.xyz")), OpenBabelNN())},
    "carbonate": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "carbonate_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "carbonate_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "carbonate_pro_-1.xyz")), OpenBabelNN())},
    "diazonium": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "diazonium_rct_1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "diazonium_ts_1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "diazonium_pro_1.xyz")), OpenBabelNN())},
    # "ester": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "ester_rct_-1.xyz")), OpenBabelNN()),
    #           "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "ester_ts_-1.xyz")), OpenBabelNN()),
    #           "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "ester_pro_-1.xyz")), OpenBabelNN())},
    "imine": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "imine_rct_0.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "imine_ts_0.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir_mols, "imine_pro_0.xyz")), OpenBabelNN())}
}

db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

mp2_tsopt = list(db.db["tasks"].find({"tags.set": "20220219_mp2_tsopt"}))
mp2_qirc = list(db.db["tasks"].find({"tags.set": "20220219_mp2_qirc"}))

grouped_data = {rxn: {"rct": dict(), "ts": dict(), "pro": dict()} for rxn in baseline_mg}

for calc in mp2_tsopt:
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

    this_rxn = None

    for rxnname, mgs in baseline_mg.items():
        mg = mgs["ts"]
        if mg.molecule.composition.alphabetical_formula == mol.composition.alphabetical_formula:
            this_rxn = rxnname
            grouped_data[rxnname]["ts"]["molecule"] = mol
            grouped_data[rxnname]["ts"]["frequencies"] = calc["output"]["frequencies"]
            grouped_data[rxnname]["ts"]["energy_svpd"] = calc["output"]["final_energy"]
            grouped_data[rxnname]["ts"]["enthalpy"] = calc["output"]["enthalpy"]
            grouped_data[rxnname]["ts"]["entropy"] = calc["output"]["entropy"]
            break
        if this_rxn is not None:
            break


for calc in mp2_qirc:
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
    calc_mg = MoleculeGraph.with_local_env_strategy(
        mol,
        OpenBabelNN()
    )

    this_rxn = None
    this_state = None

    for rxnname, mgs in baseline_mg.items():
        for state in ["rct", "pro"]:
            mg = mgs[state]
            if calc_mg.isomorphic_to(mg):
                this_rxn = rxnname
                this_state = state
                grouped_data[rxnname][state]["molecule"] = mol
                grouped_data[rxnname][state]["frequencies"] = calc["output"]["frequencies"]
                grouped_data[rxnname][state]["energy_svpd"] = calc["output"]["final_energy"]
                grouped_data[rxnname][state]["enthalpy"] = calc["output"]["enthalpy"]
                grouped_data[rxnname][state]["entropy"] = calc["output"]["entropy"]
                break
        if this_rxn is not None and this_state is not None:
            break


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

alpha_svp_tzvp = 10.390
alpha_tzvpp_qzvpp = 7.880

beta_svp_tzvp = 2.400
beta_tzvpp_qzvpp = 2.970

def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)

base_dir = "/Users/ewcss/data/ssbt/20220228_mp2_sp/block_2022-03-01-01-01-30-654215"

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

    this_rxn = None
    this_state = None

    if "perturb" in name:
        for rxnname, mgs in baseline_mg.items():
            for state in ["rct", "pro"]:
                mg = mgs[state]
                if calc_mg.isomorphic_to(mg):
                    this_rxn = rxnname
                    this_state = state
                    break
            if this_rxn is not None and this_state is not None:
                break
    else:
        for rxnname, mgs in baseline_mg.items():
            mg = mgs["ts"]
            if mg.molecule.composition.alphabetical_formula == mol.composition.alphabetical_formula:
                this_rxn = rxnname
                this_state = "ts"
                break

    with zopen(os.path.join(base_dir, direct, filename), mode="rt", encoding="ISO-8859-1") as f:
        text = f.read()

        if is_mp2:
            thishf = hf_mp2.search(text)
            if thishf is not None:
                grouped_data[this_rxn][this_state]["hf_" + basis] = float(thishf.group(1))

            thismp2_corr = mp2_mp2.search(text)
            if thismp2_corr is not None:
                grouped_data[this_rxn][this_state]["mp2_corr_" + basis] = float(thismp2_corr.group(1))

            thismp2 = mp2_total.search(text)
            if thismp2 is not None:
                grouped_data[this_rxn][this_state]["mp2_total_" + basis] = float(thismp2.group(1))

        else:
            thishf = hf.search(text)
            if thishf is not None:
                grouped_data[this_rxn][this_state]["hf_" + basis] = float(thishf.group(1))

            thismp2 = mp2.search(text)
            if thismp2 is not None:
                grouped_data[this_rxn][this_state]["ccsdt_mp2_" + basis] = float(thismp2.group(1))

            thisccsd_corr = ccsd_corr.search(text)
            if thisccsd_corr is not None:
                grouped_data[this_rxn][this_state]["ccsd_corr_" + basis] = float(thisccsd_corr.group(1))

            thisccsd_total = ccsd_total.search(text)
            if thisccsd_total is not None:
                grouped_data[this_rxn][this_state]["ccsd_total_" + basis] = float(thisccsd_total.group(1))

            thisccsdt_corr = ccsdt_corr.search(text)
            if thisccsdt_corr is not None:
                grouped_data[this_rxn][this_state]["ccsdt_corr_" + basis] = float(thisccsdt_corr.group(1))

            thisccsdt_total = ccsdt_total.search(text)
            if thisccsdt_total is not None:
                grouped_data[this_rxn][this_state]["ccsdt_total_" + basis] = float(thisccsdt_total.group(1))

full_corrected_data = dict()

for groupname, groupdata in grouped_data.items():
    for statename, statedata in groupdata.items():
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

        diff = statedata["ccsdt_total_def2-TZVP"] - statedata["ccsdt_mp2_def2-TZVP"]

        grouped_data[groupname][statename]["total_corrected"] = total_e + diff

dumpfn(grouped_data, "mp2opt_svpd_reference_data.json")