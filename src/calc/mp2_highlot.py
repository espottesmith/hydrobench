import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.qchem.outputs import QCOutput

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import (SinglePointFW,
                                          OptimizeFW,
                                          FrequencyFlatteningTransitionStateFW,
                                          TransitionStateFW,
                                          FrequencyFW)


from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

# finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220226_mp2_tzvppd"]}})]
#
# for calc in db.db["tasks"].find({"tags.set": "20220219_mp2_tsopt"}):
#     if "ester" in calc["task_label"]:
#         continue
#     mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
#
#     name = calc["task_label"].replace("SVP", "TZVPPD")
#
#     if name in finished:
#         continue
#
#     params = {"dft_rung": 4,
#               "basis_set": "def2-tzvppd",
#               "overwrite_inputs": {"rem": {"method": "mp2",
#                                            "thresh": "14",
#                                            "scf_algorithm": "diis"}}}
#
#     fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
#                                               name=name,
#                                               qchem_input_params=params,
#                                               db_file=">>db_file<<"
#                                               )
#     wf = Workflow([fw], name=name)
#     wf = add_tags(wf, {"class": "ssbt", "set": "20220226_mp2_tzvppd"})
#
#     print(wf)
#     lp.add_wf(wf)

params = {"dft_rung": 4,
          "basis_set": "def2-tzvppd",
          "overwrite_inputs": {"rem": {"method": "mp2",
                                       "thresh": "14",
                                       "scf_algorithm": "diis"}}}

# diazonium = QCOutput("/Users/ewcss/data/ssbt/20220226_mp2_tzvppd/block_2022-02-26-21-58-49-515794/launcher_2022-02-26-23-21-26-653510/ts.out")
# freq_fw = FrequencyFW(diazonium.data["molecule_from_last_geometry"],
#                       name="diazonium_ts_1_MP2_def2-TZVPPD_vacuum_freq",
#                       qchem_input_params=copy.deepcopy(params),
#                       db_file=">>db_file<<"
#                       )
# freq_wf = Workflow([freq_fw], name="diazonium_ts_1_MP2_def2-TZVPPD_vacuum_freq")
# freq_wf = add_tags(freq_wf, {"class": "ssbt", "set": "20220226_mp2_tzvppd_gaps"})
# print(freq_wf)
# lp.add_wf(freq_wf)

amide = QCOutput("/Users/ewcss/data/ssbt/20220226_mp2_tzvppd/block_2022-02-26-21-58-49-515794/launcher_2022-02-26-23-21-22-913909/ts.out")
ts_fw = FrequencyFlatteningTransitionStateFW(amide.data["molecule_from_last_geometry"],
                                             name="amide_2_2_ts_-1_MP2_def2-TZVPPD_vacuum_tsopt",
                                             qchem_input_params=copy.deepcopy(params),
                                             freq_before_opt=False,
                                             db_file=">>db_file<<"
                                             )
ts_wf = Workflow([ts_fw], name="amide_2_2_ts_-1_MP2_def2-TZVPPD_vacuum_tsopt")
ts_wf = add_tags(ts_wf, {"class": "ssbt", "set": "20220226_mp2_tzvppd_gaps"})
print(ts_wf)
lp.add_wf(ts_wf)