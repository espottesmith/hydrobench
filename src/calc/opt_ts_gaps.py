import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW, FrequencyFlatteningTransitionStateFW
from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220226_opt_gaps_tzvppd"]}})]

# choices = [("imine", "M06-HF"),
#            ("imine", "M11")]

# imine_ts = Molecule.from_file("/Users/ewcss/data/ssbt/for_sp/imine_ts_0.xyz")

# params_1 = {"dft_rung": 4,
#             "basis_set": "def2-svpd",
#             "overwrite_inputs": {"rem": {
#                 "method": "M06-HF",
#                 "thresh": "14",
#                 "scf_algorithm": "diis",
#                 "geom_opt_dmax": 50}}
#             }
#
# fw_1 = FrequencyFlatteningTransitionStateFW(molecule=copy.deepcopy(imine_ts),
#                                           name="imine_ts_0_M06-HF_def2-SVPD_vacuum_tsopt",
#                                           qchem_input_params=params_1,
#                                           db_file=">>db_file<<"
#                                           )
# wf_1 = Workflow([fw_1], name="imine_ts_0_M06-HF_def2-SVPD_vacuum_tsopt")
# wf_1 = add_tags(wf_1, {"class": "ssbt", "set": "20220226_opt_gaps"})
#
# print(wf_1)
# lp.add_wf(wf_1)
#
# params_2 = {"dft_rung": 4,
#             "basis_set": "def2-svpd",
#             "overwrite_inputs": {"rem": {
#                 "method": "M11",
#                 "thresh": "14",
#                 "scf_algorithm": "diis",
#                 "geom_opt_dmax": 50}}
#             }
#
# fw_2 = FrequencyFlatteningTransitionStateFW(molecule=copy.deepcopy(imine_ts),
#                                           name="imine_ts_0_M11_def2-SVPD_vacuum_tsopt",
#                                           qchem_input_params=params_2,
#                                           db_file=">>db_file<<"
#                                           )
# wf_2 = Workflow([fw_2], name="imine_ts_0_M11_def2-SVPD_vacuum_tsopt")
# wf_2 = add_tags(wf_2, {"class": "ssbt", "set": "20220226_opt_gaps"})
#
# print(wf_2)
# lp.add_wf(wf_2)

# for calc in db.db["tasks"].find({"tags.set": "20220130_opt_benchmark_constopts"}):
#     match = False
#     for c in choices:
#         if all([x in calc["task_label"] for x in c]):
#             match = True
#             break
#
#     if not match:
#         continue
#
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
#
#     name = calc["task_label"].replace("constopt", "tsopt")
#
#     if name in finished:
#         continue
#
#     params = {"dft_rung": 4,
#               "basis_set": "def2-svpd",
#               "overwrite_inputs": {"rem": {"thresh": "14",
#                                            "scf_algorithm": "diis",
#                                            "geom_opt_dmax": 50}}}
#
#     params["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
#     if calc["orig"]["rem"].get("dft_d"):
#         params["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]
#
#     if "IEF-PCM" in name:
#         continue
#         # params["pcm_dielectric"] = 78.39
#         # params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
#
#     fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
#                                               name=name,
#                                               qchem_input_params=params,
#                                               db_file=">>db_file<<"
#                                               )
#     wf = Workflow([fw], name=name)
#     wf = add_tags(wf, {"class": "ssbt", "set": "20220226_opt_gaps"})
#
#     print(wf)
#     lp.add_wf(wf)

# for calc in db.db["tasks"].find({"tags.set": "20220226_opt_gaps", "state": "successful"}):
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
#     mode = calc["output"]["frequency_modes"][0]
#
#     name = calc["task_label"].replace("tsopt", "")
#
#     if any([name in x for x in finished]):
#         continue
#
#     params_qirc = {"dft_rung": 4,
#                   "basis_set": "def2-svpd",
#                   "overwrite_inputs": {"rem": {"thresh": "14",
#                                                "scf_algorithm": "diis"}}}
#
#     params_qirc["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
#     if calc["orig"]["rem"].get("dft_d"):
#         params_qirc["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]
#
#     wf_2 = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
#                                         qchem_input_params=params_qirc,
#                                         name=name + "qirc")
#     wf_2 = add_tags(wf_2, {"class": "ssbt", "set": "20220226_opt_gaps_qirc"})
#
#     print(wf_2)
#     lp.add_wf(wf_2)

# for calc in db.db["tasks"].find({"tags.set": "20220226_opt_gaps"}):
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
#
#     name = calc["task_label"].replace("tsopt", "")
#
#     if any([name in x for x in finished]):
#         continue
#
#     params_tzvppd = {"dft_rung": 4,
#                   "basis_set": "def2-tzvppd",
#                   "overwrite_inputs": {"rem": {"thresh": "14",
#                                                "scf_algorithm": "diis"}}}
#
#     params_tzvppd["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
#     if calc["orig"]["rem"].get("dft_d"):
#         params_tzvppd["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]
#
#     fw = FrequencyFlatteningTransitionStateFW(molecule=copy.deepcopy(mol),
#                                               name=name + "tsopt_tzvppd",
#                                               qchem_input_params=params_tzvppd,
#                                               db_file=">>db_file<<"
#                                               )
#     wf = Workflow([fw], name=name + "tsopt_tzvppd")
#     wf = add_tags(wf, {"class": "ssbt", "set": "20220226_opt_gaps_tzvppd"})
#
#     print(wf)
#     lp.add_wf(wf)

for calc in db.db["tasks"].find({"tags.set": "20220226_opt_gaps_tzvppd", "state": "successful"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
    mode = calc["output"]["frequency_modes"][0]

    name = calc["task_label"].replace("tsopt", "")

    if any([name in x for x in finished]):
        continue

    params_qirc = {"dft_rung": 4,
                  "basis_set": "def2-tzvppd",
                  "overwrite_inputs": {"rem": {"thresh": "14",
                                               "scf_algorithm": "diis"}}}

    params_qirc["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
    if calc["orig"]["rem"].get("dft_d"):
        params_qirc["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]

    wf_2 = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
                                        qchem_input_params=params_qirc,
                                        name=name + "qirc")
    wf_2 = add_tags(wf_2, {"class": "ssbt", "set": "20220226_opt_gaps_tzvppd_qirc"})

    print(wf_2)
    lp.add_wf(wf_2)