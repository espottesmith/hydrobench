import copy
import os

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import FrequencyFlatteningTransitionStateFW
from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts
from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("launchpad.yaml")
db = QChemCalcDb.from_db_file("db.json")

base_dir = "../data/molecules/for_sp"

# Things we're not interested in, and things that already ran
exclude = ["amide_1_1", "amide_1_2", "amide_2_2", "carbonate", "diazonium", "imine"]

exclude_finished = ["acetic"]

for file in os.listdir(base_dir):
    if not file.endswith(".xyz"):
        continue
    elif any([x in file for x in exclude]):
        continue
    elif "ts" not in file:
        continue

    mol = Molecule.from_file(os.path.join(base_dir, file))
    mol.set_charge_and_spin(int(file.split(".")[0].split("_")[-1]))

    params = {"basis_set": "def2-svpd",
              "overwrite_inputs": {"rem": {"method": "wB97M-V",
                                           "scf_algorithm": "diis",
                                           "thresh": 14}
                                   }
              }

    name = file.split(".")[0] + "fftsopt"

    fw = FrequencyFlatteningTransitionStateFW(molecule=mol,
                                              name=name,
                                              qchem_input_params=params,
                                              db_file=">>db_file<<"
                                              )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20220305_mv_fftsopt"})

    print(wf)
    lp.add_wf(wf)

for calc in db.db["tasks"].find({"tags.set": "20220305_mv_fftsopt"}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])
    mode = calc["output"]["frequency_modes"][0]
    name = calc["task_label"]

    if not any([x in name for x in exclude_finished]):
        continue

    params = {"basis_set": "def2-svpd",
              "overwrite_inputs": {"rem": {"method": "wB97M-V",
                                           "scf_algorithm": "diis",
                                           "thresh": 14}
                                   }
              }

    wf = get_wf_reaction_path_with_ts(copy.deepcopy(mol), mode, name, scale=0.6,
                                      qchem_input_params=params,
                                      name=name + "qirc")
    wf = add_tags(wf, {"class": "ssbt", "set": "20220308_mv_qirc"})

    print(wf)
    lp.add_wf(wf)