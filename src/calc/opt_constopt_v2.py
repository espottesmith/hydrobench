import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts
from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW, FrequencyFlatteningOptimizeFW, FrequencyFlatteningTransitionStateFW

from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

dft_methods = ["PBE", ("PBE", "D3(BJ)"), "BLYP", ("BLYP", "D3(BJ)"), "B97-D", "B97-D3", "mPW91", ("mPW91", "D3(BJ)"), "VV10", "rVV10",
               "M06-L", ("M06-L", "D3(0)"), "SCAN", ("SCAN", "D3(BJ)"), "TPSS", ("TPSS", "D3(BJ)"), "MN12-L", ("MN12-L", "D3(BJ)"), "B97M-rV",
               "PBE0", ("PBE0", "D3(BJ)"), "LRC-wPBE", ("LRC-wPBE", "D3(BJ)"), "LRC-wPBEh", ("LRC-wPBEh", "D3(BJ)"), "B3LYP", ("B3LYP", "D3(BJ)"), "CAM-B3LYP", ("CAM-B3LYP", "D3(0)"), "rCAM-B3LYP", ("rCAM-B3LYP", "D3(0)"), "mPW1PW91", ("mPW1PW91", "D3(BJ)"), "HSE-HJS", ("HSE-HJS", "D3(BJ)"), "wB97X", "wB97XD", "wB97XD3", "wB97XV",
               "M06-2X", ("M06-2X", "D3(0)"), "wM06-D3", "M06-SX", ("M06-SX", "D3(BJ)"), "M06-HF", ("M06-HF", "D3(0)"), "M08-SO", ("M08-SO", "D3(0)"), "M11", ("M11", "D3(0)"), "revM11", ("revM11", "D3(0)"), "MN15", ("MN15", "D3(0)"), "BMK", ("BMK", "D3(BJ)"), "TPSSh", ("TPSSh", "D3(BJ)"), "SCAN0", ("SCAN0", "D3(BJ)"), "mPWB1K", ("mPWB1K", "D3(BJ)"), "wB97M-V"]

# Eventually reopt in def2-TZVPPD
dft_basis = "def2-SVPD"

# constraints = {
#     "aceticanhydride_ts_0": ["6 7 0.973", "4 7 1.876", "1 2 1.405", "1 6 2.866"],
#     "amide_2_1_ts_-1": ["2 13 1.895"],
#     "amide_2_2_ts_-1": ["2 9 1.834", "4 8 1.005", "8 16 1.712", "9 12 1.208", "12 16 1.312"],
#     "borohydride_ts_-1": ["4 7 1.159", "5 7 1.306", "1 5 2.196"],
#     "carbonate_ts_-1": ["1 5 1.738", "1 7 2.110"],
#     "diazonium_ts_1": ["4 12 2.255", "4 14 2.438"],
#     "basic_epoxide_1_ts_-1": ["1 3 1.655", "1 8 2.299", "2 3 1.409"],
#     "basic_epoxide_2_ts_-1": ["1 3 1.719", "1 8 2.290", "2 3 1.400"],
#     "epoxide1_ts_1": ["1 3 1.881", "1 10 2.315", "2 3 1.457"],
#     "epoxide2_ts_1": ["1 3 1.807", "1 10 2.102", "2 3 1.484"],
#     "ester_ts_-1": ["1 7 2.628", "4 7 2.110", "3 4 1.110"],
#     "furan1_ts_1": ["1 3 2.15", "3 17 1.95", "17 19 1.002", "19 20 1.610"],
#     "furan2_ts_1": ["3 17 2.111", "17 19 0.992", "19 20 1.671"],
#     "furan3_ts_1": ["5 21 1.621", "20 21 1.104", "17 22 1.505", "3 17 1.435"],
#     "imine_ts_0": ["2 8 1.062", "3 4 2.273", "4 8 1.705"],
#     "iminium_ts_1": ["2 8 1.461", "4 8 1.158", "3 4 1.498"],
#     "lactone_ts_-1": ["1 3 1.473", "1 9 1.625"],
# constraints = {"cl2co_ts_0": ["1 5 1.592", "5 6 1.195", "2 6 1.318"]}

base_dir = "/Users/ewcss/data/ssbt/for_sp"

# finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20230403_opt_benchmark_constopt"]}})]

# for mol_name, const in constraints.items():
#     for mol_file in os.listdir(base_dir):
#         if mol_name in mol_file:
#             mol = Molecule.from_file(os.path.join(base_dir, mol_file))
#             mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

#             for method in dft_methods:
#                 params = {"basis_set": dft_basis,
#                           "qchem_version": 5,
#                           "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
#                                                        "thresh": 14},
#                                                "opt": {"CONSTRAINT": ["stre " + x for x in const]}}}
#                 if isinstance(method, tuple):
#                     params["overwrite_inputs"]["rem"]["method"] = method[0]
#                     if method[1] == "D3(BJ)":
#                         params["overwrite_inputs"]["rem"]["dft_d"] = "D3_BJ"
#                     elif method[1] == "D3(0)":
#                         params["overwrite_inputs"]["rem"]["dft_d"] = "D3_ZERO"

#                     base_name = mol_file.split(".")[0] + "_" + method[0] + "-" + method[1] + "_" + dft_basis
#                 else:
#                     params["overwrite_inputs"]["rem"]["method"] = method
#                     base_name = mol_file.split(".")[0] + "_" + method + "_" + dft_basis

#                 name_vac = base_name + "_" + "vacuum"
#                 if name_vac not in finished:
#                     fw_vac = OptimizeFW(copy.deepcopy(mol),
#                                         name=name_vac + "_constopt",
#                                         qchem_input_params=copy.deepcopy(params),
#                                         db_file=">>db_file<<")
#                     wf_vac = Workflow([fw_vac], name=name_vac + "_constopt")
#                     wf_vac = add_tags(wf_vac, {"class": "ssbt",
#                                                "set": "20230403_opt_benchmark_constopt"})
#                     print(wf_vac)
#                     lp.add_wf(wf_vac)


# finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20230407_opt_benchmark_tsopt"]}})]

# for calc in db.db["tasks"].find({"tags.set": "20230403_opt_benchmark_constopt"}):
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

#     name = calc["task_label"].replace("constopt", "tsopt")

#     if name in finished:
#         continue

#     params = {"dft_rung": 4,
#               "basis_set": "def2-svpd",
#               "overwrite_inputs": {"rem": {"thresh": "14",
#                                            "scf_algorithm": "diis"}}}

#     params["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
#     if calc["orig"]["rem"].get("dft_d"):
#         params["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]

#     fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
#                                               name=name,
#                                               qchem_input_params=params,
#                                               db_file=">>db_file<<"
#                                               )
#     wf = Workflow([fw], name=name)
#     wf = add_tags(wf, {"class": "ssbt", "set": "20230407_opt_benchmark_tsopt"})

#     print(wf)
#     lp.add_wf(wf)


finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20230408_opt_benchmark_qirc"]}})]
for calc in db.db["tasks"].find({"tags.set": "20230407_opt_benchmark_tsopt"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
    mode = calc["output"]["frequency_modes"][0]

    name = calc["task_label"].replace("tsopt", "")

    if any([name in x for x in finished]):
        continue

    params_qirc = {"qchem_version": 5,
                   "basis_set": "def2-svpd",
                   "overwrite_inputs": {"rem": {"thresh": "14",
                                                "scf_algorithm": "diis"}}}

    params_qirc["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
    if calc["orig"]["rem"].get("dft_d"):
        params_qirc["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]

    if "IEF-PCM" in name:
        continue

    wf_2 = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
                                        qchem_input_params=params_qirc,
                                        name=name + "qirc")
    wf_2 = add_tags(wf_2, {"class": "ssbt", "set": "20230408_opt_benchmark_qirc"})

    print(wf_2)
    lp.add_wf(wf_2)

# for mol_file in os.listdir(base_dir):
#     if not mol_file.endswith(".xyz"):
#         continue

#     if "_ts_" not in mol_file:
#         mol = Molecule.from_file(os.path.join(base_dir, mol_file))
#         mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

#         for method in dft_methods:
#             params = {"basis_set": dft_basis,
#                       "qchem_version": 5,
#                       "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
#                                                    "thresh": 14}
#                                           }
#                     }
#             if isinstance(method, tuple):
#                 params["overwrite_inputs"]["rem"]["method"] = method[0]
#                 if method[1] == "D3(BJ)":
#                     params["overwrite_inputs"]["rem"]["dft_d"] = "D3_BJ"
#                 elif method[1] == "D3(0)":
#                     params["overwrite_inputs"]["rem"]["dft_d"] = "D3_ZERO"

#                 base_name = mol_file.split(".")[0] + "_" + method[0] + "-" + method[1] + "_" + dft_basis
#             else:
#                 params["overwrite_inputs"]["rem"]["method"] = method
#                 base_name = mol_file.split(".")[0] + "_" + method + "_" + dft_basis

#             name_vac = base_name + "_" + "vacuum"
#             # if name_vac not in finished:
#             fw_vac = FrequencyFlatteningOptimizeFW(copy.deepcopy(mol),
#                                 name=name_vac + "_directopt",
#                                 qchem_input_params=copy.deepcopy(params),
#                                 db_file=">>db_file<<")
#             wf_vac = Workflow([fw_vac], name=name_vac + "_directopt")
#             wf_vac = add_tags(wf_vac, {"class": "ssbt",
#                                         "set": "20230408_opt_benchmark_directopt"})
#             print(wf_vac)
#             lp.add_wf(wf_vac)