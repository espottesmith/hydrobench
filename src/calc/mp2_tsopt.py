from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import FrequencyFlatteningTransitionStateFW
from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("launchpad.yaml")
db = QChemCalcDb.from_db_file("db.json")

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220219_mp2_tsopt"]}})]

for calc in db.db["tasks"].find({"tags.set": "20220218_mp2_constopt"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

    name = calc["task_label"].replace("constopt", "tsopt")

    if name in finished:
        continue

    params = {"dft_rung": 4,
              "basis_set": "def2-svpd",
              "overwrite_inputs": {"rem": {"method": "mp2",
                                           "thresh": "14",
                                           "scf_algorithm": "diis"}}}

    fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
                                              name=name,
                                              qchem_input_params=params,
                                              db_file=">>db_file<<"
                                              )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20220219_mp2_tsopt"})

    print(wf)
    lp.add_wf(wf)