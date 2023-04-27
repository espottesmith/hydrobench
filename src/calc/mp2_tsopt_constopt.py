import copy
import os

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import OptimizeFW
from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

basis = "def2-SVP"

constraints = {
    "aceticanhydride_ts_0": ["6 7 0.973", "4 7 1.876", "1 2 1.405", "1 6 2.866"],
    "amide_2_1_ts_-1": ["2 13 1.895"],
    "amide_2_2_ts_-1": ["2 9 1.834", "4 8 1.005", "8 16 1.712", "9 12 1.208", "12 16 1.312"],
    "borohydride_ts_-1": ["4 7 1.159", "5 7 1.306", "1 5 2.196"],
    "carbonate_ts_-1": ["1 5 1.738", "1 7 2.110"],
    "diazonium_ts_1": ["4 12 2.255", "4 14 2.438"],
    "basic_epoxide_1_ts_-1": ["1 3 1.655", "1 8 2.299", "2 3 1.409"],
    "basic_epoxide_2_ts_-1": ["1 3 1.719", "1 8 2.290", "2 3 1.400"],
    "epoxide1_ts_1": ["1 3 1.881", "1 10 2.315", "2 3 1.457"],
    "epoxide2_ts_1": ["1 3 1.807", "1 10 2.102", "2 3 1.484"],
    "ester_ts_-1": ["1 7 2.628", "4 7 2.110", "3 4 1.110"],
    "furan1_ts_1": ["1 3 2.15", "3 17 1.95", "17 19 1.002", "19 20 1.610"],
    "furan2_ts_1": ["3 17 2.111", "17 19 0.992", "19 20 1.671"],
    "furan3_ts_1": ["5 21 1.621", "20 21 1.104", "17 22 1.505", "3 17 1.435"],
    "imine_ts_0": ["2 8 1.062", "3 4 2.273", "4 8 1.705"],
    "iminium_ts_1": ["2 8 1.461", "4 8 1.158", "3 4 1.498"],
    "lactone_ts_-1": ["1 3 1.473", "1 9 1.625"],
    "cl2co_ts_0": ["1 5 1.592", "5 6 1.195", "2 6 1.318"]}

base_dir = "../data/molecules/for_sp"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220130_opt_benchmark_constopts"]}})]

for mol_name, const in constraints.items():
    for mol_file in os.listdir(base_dir):
        if mol_name in mol_file:
            mol = Molecule.from_file(os.path.join(base_dir, mol_file))
            mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

            params = {"basis_set": basis,
                      "overwrite_inputs": {"rem": {
                                                "method": "mp2",
                                                "mem_total": 760000,
                                                "cc_memory": 600000,
                                            },
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
                                           "set": "20230411_mp2_constopt"})
                print(wf_vac)
                lp.add_wf(wf_vac)