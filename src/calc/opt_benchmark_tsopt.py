import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW, FrequencyFlatteningTransitionStateFW

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220130_opt_benchmark_tsopt", "20220203_opt_benchmark_tsopt"]}})]

for calc in db.db["tasks"].find({"tags.set": "20220130_opt_benchmark_constopts"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

    name = calc["task_label"].replace("constopt", "tsopt")

    if name in finished:
        continue

    params = {"dft_rung": 4,
              "basis_set": "def2-svpd",
              "overwrite_inputs": {"rem": {"thresh": "14",
                                           "scf_algorithm": "diis"}}}

    params["overwrite_inputs"]["rem"]["method"] = calc["orig"]["rem"]["method"]
    if calc["orig"]["rem"].get("dft_d"):
        params["overwrite_inputs"]["rem"]["dft_d"] = calc["orig"]["rem"]["dft_d"]

    if "IEF-PCM" in name:
        continue
        # params["pcm_dielectric"] = 78.39
        # params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}

    fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
                                              name=name,
                                              qchem_input_params=params,
                                              db_file=">>db_file<<"
                                              )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20220203_opt_benchmark_tsopt"})

    print(wf)
    lp.add_wf(wf)