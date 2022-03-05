import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW, FrequencyFlatteningTransitionStateFW
from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts

from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220207_opt_benchmark_tzvppd_sp", "20220227_opt_benchmark_tzvppd_sp"]}})]

# for calc in db.db["tasks"].find({"tags.set": "20220203_opt_benchmark_tsopt"}):
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
#     if "IEF-PCM" in name:
#         # params_qirc["pcm_dielectric"] = 78.39
#         # params_qirc["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
#         continue
#
#     params_tzvppd = copy.deepcopy(params_qirc)
#     params_tzvppd["basis_set"] = "def2-tzvppd"
#
#     fw = FrequencyFlatteningTransitionStateFW(molecule=copy.deepcopy(mol),
#                                               name=name + "tsopt_tzvppd",
#                                               qchem_input_params=params_tzvppd,
#                                               db_file=">>db_file<<"
#                                               )
#     wf = Workflow([fw], name=name + "tsopt_tzvppd")
#     wf = add_tags(wf, {"class": "ssbt", "set": "20220204_opt_benchmark_tzvppd"})
#
#     print(wf)
#     lp.add_wf(wf)
#
#     wf_2 = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
#                                         qchem_input_params=params_qirc,
#                                         name=name + "qirc")
#     wf_2 = add_tags(wf_2, {"class": "ssbt", "set": "20220204_opt_benchmark_qirc"})
#
#     print(wf_2)
#     lp.add_wf(wf_2)


# for calc in db.db["tasks"].find({"tags.set": "20220204_opt_benchmark_tzvppd"}):
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
#     mode = calc["output"]["frequency_modes"][0]
#
#     name = calc["task_label"].replace("tsopt", "")
#
#     if any([name in x for x in finished]):
#         continue
#
#     params = {"dft_rung": 4,
#               "basis_set": "def2-tzvppd",
#               "overwrite_inputs": {"rem": {"thresh": "14",
#                                            "scf_algorithm": "diis"}}}
#
#     params["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
#     if calc["orig"]["rem"].get("dft_d"):
#         params["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]
#
#     if "IEF-PCM" in name:
#         # params_qirc["pcm_dielectric"] = 78.39
#         # params_qirc["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
#         continue
#
#     wf = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
#                                         qchem_input_params=params,
#                                         name=name + "qirc")
#     wf = add_tags(wf, {"class": "ssbt", "set": "20220207_opt_benchmark_tzvppd_qirc"})
#
#     print(wf)
#     lp.add_wf(wf)

for calc in db.db["tasks"].find({"tags.set": {"$in": ["20220226_opt_gaps", "20220226_opt_gaps_qirc"]}}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

    name = calc["task_label"] + "SP"

    # if any([name in x for x in finished]):
    #     continue

    params = {"dft_rung": 4,
              "basis_set": "def2-tzvppd",
              "overwrite_inputs": {"rem": {"thresh": "14",
                                           "scf_algorithm": "diis"}}}

    params["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
    if calc["orig"]["rem"].get("dft_d"):
        params["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]

    if "IEF-PCM" in name:
        continue

    fw = SinglePointFW(molecule=mol,
                       name=name,
                       qchem_input_params=params,
                       db_file=">>db_file<<"
                       )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20220227_opt_benchmark_tzvppd_sp"})
    print(wf)
    lp.add_wf(wf)