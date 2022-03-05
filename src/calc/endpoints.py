import copy

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

params = {"dft_rung": 4,
          "basis_set": "def2-svpd",
          "overwrite_inputs": {"rem": {"method": "scan",
                                       "thresh": "14",
                                       "scf_algorithm": "diis",
                                       "geom_opt_dmax": 100}}}

ts_to_qirc = db.db["tasks"].find({"tags.class": "ssbt", "special_run_type": "ts_frequency_flattener", "state": "successful", "output.frequencies.0": {"$lt": 0}, "output.frequencies.1": {"$gt": 0}})

already_done = [e["task_label"] for e in db.db["tasks"].find({"tags.set": "20211004_qirc"})]

for ts in ts_to_qirc:
    mol = Molecule.from_dict(ts["output"]["optimized_molecule"])
    mol.to("xyz", "/Users/ewcss/data/ssbt/20220118_ts_opted/{}.xyz".format(ts["task_label"]))
    mode = ts["output"]["frequency_modes"][0]
    name = ts["task_label"].split()[0] + "_qirc"

    match = False
    for a in already_done:
        if name in a:
            match = True
            break

    if match:
        continue

    wf = get_wf_reaction_path_with_ts(mol, mode, name, scale=0.6, qchem_input_params=copy.deepcopy(params),
                                      name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20211004_qirc"})

    print(wf)
    lp.add_wf(wf)