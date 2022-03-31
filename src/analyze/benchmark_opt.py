import copy
import os
import statistics

import numpy as np

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.molecule_matcher import MoleculeMatcher, IsomorphismMolAtomMapper

from atomate.qchem.database import QChemCalcDb


def perturb(mol, mode, scale=0.6):
    mol_copy = copy.deepcopy(mol)
    for ii in range(len(mol)):
        vec = np.array(mode[ii])
        mol_copy.translate_sites(indices=[ii], vector=vec * scale)
    return mol_copy


def get_free_energy(energy, enthalpy, entropy, temp=298.15):
    return energy * 27.2114 + enthalpy * 0.043363 - temp * entropy * 0.000043363


# Associate TS with reactants and products (separately with TZVPPD)
# First, just make sure that all of the calculations wind up with the same reactants and products
# Then, look at active bonds at the TS and look for differences
# And maybe look at difference between barriers with SCAN geometry and with functional-optimized geometry?

db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

base_dir = "/Users/ewcss/data/ssbt/for_sp"

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

baseline_mg = {
    "amide": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_pro_-1.xyz")), OpenBabelNN())},
    "carbonate": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_pro_-1.xyz")), OpenBabelNN())},
    "diazonium": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_rct_1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_ts_1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_pro_1.xyz")), OpenBabelNN())},
    # "ester": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_rct_-1.xyz")), OpenBabelNN()),
    #           "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_ts_-1.xyz")), OpenBabelNN()),
    #           "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_pro_-1.xyz")), OpenBabelNN())},
    "imine": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_rct_0.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_ts_0.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_pro_0.xyz")), OpenBabelNN())}
}

svpd_ts = list(db.db["tasks"].find({"tags.set": {"$in": ["20220203_opt_benchmark_tsopt", "20220226_opt_gaps"]}}))
svpd_end = list(db.db["tasks"].find({"tags.set": {"$in": ["20220204_opt_benchmark_qirc", "20220226_opt_gaps_qirc"]}}))
tzvppd_sp = list(db.db["tasks"].find({"tags.set": {"$in": ["20220207_opt_benchmark_tzvppd_sp", "20220227_opt_benchmark_tzvppd_sp"]}}))

# tzvppd_ts = list(db.db["tasks"].find({"tags.set": {"$in": ["20220204_opt_benchmark_tzvppd", "20220226_opt_gaps_tzvppd"]}}))
# tzvppd_end = list(db.db["tasks"].find({"tags.set": {"$in": ["20220207_opt_benchmark_tzvppd_qirc", "20220226_opt_gaps_tzvppd_qirc"]}}))

data_sets = {n: dict() for n in baseline_mg}

no_endpoints = list()
different_endpoints = list()

for ts_calc in svpd_ts:
    if "PCM" in ts_calc["task_label"]:
        continue

    elif "ester" in ts_calc["task_label"]:
        continue

    ts_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(ts_calc["output"]["optimized_molecule"]), OpenBabelNN())
    mode = ts_calc["output"]["frequency_modes"][0]

    name = copy.deepcopy(ts_calc["task_label"])
    for a, b in replacements.items():
        name = name.replace(a, b)

    contents = name.split("_")
    method = contents[3]

    for_mol = perturb(ts_mg.molecule, mode, scale=0.6)
    rev_mol = perturb(ts_mg.molecule, mode, scale=-0.6)

    for_calc = None
    for_mg = None
    rev_calc = None
    rev_mg = None

    for end_calc in svpd_end:
        if ts_calc["task_label"].replace("tsopt", "") not in end_calc["task_label"]:
            continue
        end_mol = Molecule.from_dict(end_calc["orig"]["molecule"])
        if end_mol == for_mol and for_calc is None:
            for_calc = end_calc
            for_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(end_calc["output"]["optimized_molecule"]),
                                                           OpenBabelNN())
        elif end_mol == rev_mol and rev_calc is None:
            rev_calc = end_calc
            rev_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(end_calc["output"]["optimized_molecule"]),
                                                           OpenBabelNN())

        if for_calc is not None and rev_calc is not None:
            break

    if for_calc is None or rev_calc is None:
        print("Endpoints could not be found for", ts_calc["task_id"], ts_calc["task_label"])
        no_endpoints.append(ts_calc["task_id"])
        continue

    this_data = {"ts": ts_calc, "rct": None, "pro": None}
    problem = False

    rct_mg = None
    pro_mg = None
    this_rxn = None
    for n, data in baseline_mg.items():
        if n in ts_calc["task_label"]:
            this_rxn = n
            if for_mg.isomorphic_to(data["pro"]) and rev_mg.isomorphic_to(data["rct"]):
                this_data["rct"] = rev_calc
                this_data["pro"] = for_calc
                rct_mg = rev_mg
                pro_mg = for_mg
            elif for_mg.isomorphic_to(data["rct"]) and rev_mg.isomorphic_to(data["pro"]):
                this_data["rct"] = for_calc
                this_data["pro"] = rev_calc
                rct_mg = for_mg
                pro_mg = rev_mg
            else:

                print("Different reactants and products for", ts_calc["task_id"], ts_calc["task_label"])
                print("FORWARD", for_calc["task_label"], "REVERSE", rev_calc["task_label"])
                different_endpoints.append(ts_calc["task_id"])
            break

    if this_data["rct"] is None or this_data["pro"] is None:
        print("No match for", ts_calc["task_id"], ts_calc["task_label"])
        continue

    for sp_calc in tzvppd_sp:
        mol = Molecule.from_dict(sp_calc["input"]["initial_molecule"])
        if mol == rct_mg.molecule:
            this_data["rct_sp"] = sp_calc
        elif mol == pro_mg.molecule:
            this_data["pro_sp"] = sp_calc
        elif mol == ts_mg.molecule:
            this_data["ts_sp"] = sp_calc
        if this_data.get("ts_sp") and this_data.get("rct_sp") and this_data.get("pro_sp"):
            break

    if not this_data.get('ts_sp'):
        print(this_rxn, method, "NO TS SP", this_data["ts"]["task_id"])
    if not this_data.get('rct_sp'):
        print(this_rxn, method, "NO RCT SP", this_data["rct"]["task_id"])
    if not this_data.get('pro_sp'):
        print(this_rxn, method, "NO PRO SP", this_data["pro"]["task_id"])

    data_sets[this_rxn][method] = this_data

# exclude = [12716, 12717, 12719, 12720]
#
# data_sets_tzvppd = {n: dict() for n in baseline_mg}
# for ts_calc in tzvppd_ts:
#
#     if "PCM" in ts_calc["task_label"]:
#         continue
#
#     elif "ester" in ts_calc["task_label"]:
#         continue
#
#     ts_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(ts_calc["output"]["optimized_molecule"]), OpenBabelNN())
#     mode = ts_calc["output"]["frequency_modes"][0]
#
#     name = copy.deepcopy(ts_calc["task_label"])
#     for a, b in replacements.items():
#         name = name.replace(a, b)
#
#     contents = name.split("_")
#     method = contents[3]
#
#     for_mol = perturb(ts_mg.molecule, mode, scale=0.6)
#     rev_mol = perturb(ts_mg.molecule, mode, scale=-0.6)
#
#     for_calc = None
#     for_mg = None
#     rev_calc = None
#     rev_mg = None
#
#     for end_calc in tzvppd_end:
#         # Cannot be the correct endpoint
#         if end_calc["task_id"] < ts_calc["task_id"]:
#             continue
#         elif end_calc["task_id"] in exclude:
#             continue
#         if ts_calc["task_label"].replace("tsopt", "") not in end_calc["task_label"]:
#             continue
#         end_mol = Molecule.from_dict(end_calc["orig"]["molecule"])
#         if (end_mol == for_mol and for_calc is None) or "forwards" in end_calc["task_label"]:
#             for_calc = end_calc
#             for_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(end_calc["output"]["optimized_molecule"]),
#                                                            OpenBabelNN())
#         elif (end_mol == rev_mol and rev_calc is None) or "backwards" in end_calc["task_label"]:
#             rev_calc = end_calc
#             rev_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(end_calc["output"]["optimized_molecule"]),
#                                                            OpenBabelNN())
#
#         if for_calc is not None and rev_calc is not None:
#             break
#
#     if for_calc is None or rev_calc is None:
#         print("Endpoints could not be found for", ts_calc["task_id"], ts_calc["task_label"])
#         no_endpoints.append(ts_calc["task_id"])
#         continue
#
#     this_data = {"ts": ts_calc, "rct": None, "pro": None}
#     problem = False
#
#     rct_mg = None
#     pro_mg = None
#     this_rxn = None
#     for n, data in baseline_mg.items():
#         if n in ts_calc["task_label"]:
#             this_rxn = n
#             if for_mg.isomorphic_to(data["pro"]) and rev_mg.isomorphic_to(data["rct"]):
#                 this_data["rct"] = rev_calc
#                 this_data["pro"] = for_calc
#                 rct_mg = rev_mg
#                 pro_mg = for_mg
#             elif for_mg.isomorphic_to(data["rct"]) and rev_mg.isomorphic_to(data["pro"]):
#                 this_data["rct"] = for_calc
#                 this_data["pro"] = rev_calc
#                 rct_mg = for_mg
#                 pro_mg = rev_mg
#             else:
#                 print("Different reactants and products for", ts_calc["task_id"], ts_calc["task_label"])
#                 print("FORWARD", for_calc["task_id"], for_calc["task_label"], "REVERSE", rev_calc["task_id"], rev_calc["task_label"])
#                 different_endpoints.append(ts_calc["task_id"])
#             break
#
#     if this_data["rct"] is None or this_data["pro"] is None:
#         continue
#
#     data_sets_tzvppd[this_rxn][method] = this_data

dumpfn(data_sets, "/Users/ewcss/data/ssbt/20220207_opt_benchmark/opt_data_dump_svpd.json")
# dumpfn(data_sets_tzvppd, "/Users/ewcss/data/ssbt/20220207_opt_benchmark/opt_data_dump_tzvppd.json")