import copy
import os

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW

from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")


dft_methods = [("PBE", "D4"), ("BLYP", "D4"), ("B97", "D4"), ("M06-L", "D4"), ("SCAN", "D4"), ("TPSS", "D4"), ("PBE0", "D4"),
               ("B3LYP", "D4"), ("CAM-B3LYP", "D4"), ("mPW1PW91", "D4"), ("wB97X", "D4"), ("TPSSh", "D4")]

dft_basis = "def2-TZVPPD"

base_dir = "../data/molecules/for_sp"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20211016_sp_dft", "20211016_sp_cc"]}})]

for mol_file in os.listdir(base_dir):

    if mol_file != "aceticanhydride_realts_0.xyz":
        continue

    mol = Molecule.from_file(os.path.join(base_dir, mol_file))
    mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

    for method in dft_methods:
        params = {"basis_set": dft_basis,
                  "qchem_version": 6,
                  "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                               "thresh": 14}}}
        if isinstance(method, tuple):
            params["overwrite_inputs"]["rem"]["method"] = method[0]
            if method[1] == "D4":
                params["overwrite_inputs"]["rem"]["dft_d"] = "D4"

            base_name = mol_file.split(".")[0] + "_" + method[0] + "-" + method[1] + "_" + dft_basis
        else:
            params["overwrite_inputs"]["rem"]["method"] = method
            base_name = mol_file.split(".")[0] + "_" + method + "_" + dft_basis

        name_vac = base_name + "_" + "vacuum"
        if name_vac not in finished:
            fw_vac = SinglePointFW(copy.deepcopy(mol),
                                   name=name_vac,
                                   qchem_input_params=copy.deepcopy(params),
                                   db_file=">>db_file<<")
            wf_vac = Workflow([fw_vac], name=name_vac)
            wf_vac = add_tags(wf_vac, {"class": "ssbt",
                                   "set": "20230702_sp_d4"})
            print(wf_vac)
            lp.add_wf(wf_vac)