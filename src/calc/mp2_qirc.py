import copy

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts
from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("launchpad.yaml")
db = QChemCalcDb.from_db_file("db.json")

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220219_mp2_qirc"]}})]


for calc in db.db["tasks"].find({"tags.set": "20220219_mp2_tsopt"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
    mode = calc["output"]["frequency_modes"][0]

    name = calc["task_label"].replace("tsopt", "")

    if any([name in x for x in finished]):
        continue

    params = {"dft_rung": 4,
              "basis_set": "def2-svpd",
              "overwrite_inputs": {"rem": {"method": "mp2",
                                           "thresh": "14",
                                           "scf_algorithm": "diis"}}}

    wf = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
                                        qchem_input_params=params,
                                        name=name + "qirc")
    wf = add_tags(wf, {"class": "ssbt", "set": "20220219_mp2_qirc"})

    print(wf)
    lp.add_wf(wf)