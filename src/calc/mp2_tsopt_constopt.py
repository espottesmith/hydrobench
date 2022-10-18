import copy
import os

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW
from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("launchpad.yaml")
db = QChemCalcDb.from_db_file("db.json")

basis = "def2-SVP"

constraints = {"diazonium_ts_1": ["4 12 2.255", "4 14 2.438"],
               "ester_ts_-1": ["1 7 2.628", "4 7 2.110", "3 4 1.110"],
               "imine_ts_0": ["2 8 1.062", "3 4 2.273", "4 8 1.705"],
               "carbonate_ts_-1": ["1 5 1.738", "1 7 2.110"],
               "amide_2_2_ts_-1": ["2 9 1.834", "4 8 1.005", "8 16 1.712", "9 12 1.208", "12 16 1.312"]}

base_dir = "../data/molecules/for_sp"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220130_opt_benchmark_constopts"]}})]

for mol_name, const in constraints.items():
    for mol_file in os.listdir(base_dir):
        if mol_name in mol_file:
            mol = Molecule.from_file(os.path.join(base_dir, mol_file))
            mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

            params = {"basis_set": basis,
                      "overwrite_inputs": {"rem": {"method": "mp2",
                                                   "scf_algorithm": "diis",
                                                   "thresh": 14},
                                           "opt": {"CONSTRAINT": ["stre " + x for x in const]}}}

            base_name = mol_file.split(".")[0] + "_" + "MP2" + "_" + basis

            name_vac = base_name + "_" + "vacuum"
            if name_vac not in finished:
                fw_vac = OptimizeFW(copy.deepcopy(mol),
                                       name=name_vac + "_constopt",
                                       qchem_input_params=copy.deepcopy(params),
                                       db_file=">>db_file<<")
                wf_vac = Workflow([fw_vac], name=name_vac + "_constopt")
                wf_vac = add_tags(wf_vac, {"class": "ssbt",
                                       "set": "20220218_mp2_constopt"})
                print(wf_vac)
                lp.add_wf(wf_vac)