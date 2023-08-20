import copy
import os

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import SinglePointFW
from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

# MP2: def2-TZVP, def2-TZVPP, def2-QZVPP
# CCSD(T): def2-TZVP


dft_methods = ["PBE", ("PBE", "D3(BJ)"), "BLYP", ("BLYP", "D3(BJ)"), "B97-D", "B97-D3", "mPW91", ("mPW91", "D3(BJ)"), "VV10", "rVV10",
               "M06-L", ("M06-L", "D3(0)"), "SCAN", ("SCAN", "D3(BJ)"), "TPSS", ("TPSS", "D3(BJ)"), "MN12-L", ("MN12-L", "D3(BJ)"), "B97M-rV",
               "PBE0", ("PBE0", "D3(BJ)"), "LRC-wPBE", ("LRC-wPBE", "D3(BJ)"), "LRC-wPBEh", ("LRC-wPBEh", "D3(BJ)"), "B3LYP", ("B3LYP", "D3(BJ)"), "CAM-B3LYP", ("CAM-B3LYP", "D3(0)"), "rCAM-B3LYP", ("rCAM-B3LYP", "D3(0)"), "mPW1PW91", ("mPW1PW91", "D3(BJ)"), "HSE-HJS", ("HSE-HJS", "D3(BJ)"), "wB97X", "wB97XD", "wB97XD3", "wB97XV",
               "M06-2X", ("M06-2X", "D3(0)"), "wM06-D3", "M06-SX", ("M06-SX", "D3(BJ)"), "M06-HF", ("M06-HF", "D3(0)"), "M08-SO", ("M08-SO", "D3(0)"), "M11", ("M11", "D3(0)"), "revM11", ("revM11", "D3(0)"), "MN15", ("MN15", "D3(0)"), "BMK", ("BMK", "D3(BJ)"), "TPSSh", ("TPSSh", "D3(BJ)"), "SCAN0", ("SCAN0", "D3(BJ)"), "mPWB1K", ("mPWB1K", "D3(BJ)"), "wB97M-V"]

dft_basis = "def2-TZVPPD"

ccsdt_bases = ["def2-TZVP"]

mp2_bases = ["def2-TZVP", "def2-TZVPP", "def2-QZVPP"]

base_dir = "../data/molecules/mv"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20230702_mv_sp_dft", "20230702_mv_sp_wf"]}})]

for mol_file in os.listdir(base_dir):

    if mol_file not in [
        "aceticanhydride_ts_0.xyz",
        "aceticanhydride_rct_0.xyz",
        "aceticanhydride_pro_0.xyz",
        "lactone_ts_-1.xyz",
        "lactone_rct_-1.xyz",
        "lactone_pro_-1.xyz",
    ]:
        continue

    try:
        mol = Molecule.from_file(os.path.join(base_dir, mol_file))
        mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))
    except:
        print(mol_file)

    for method in dft_methods:
        params = {"basis_set": dft_basis,
                  "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                               "thresh": 14}}}
        if isinstance(method, tuple):
            params["overwrite_inputs"]["rem"]["method"] = method[0]
            if method[1] == "D3(BJ)":
                params["overwrite_inputs"]["rem"]["dft_d"] = "D3_BJ"
            elif method[1] == "D3(0)":
                params["overwrite_inputs"]["rem"]["dft_d"] = "D3_ZERO"

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
                                        "set": "20230702_mv_sp_dft"})
            print(wf_vac)
            lp.add_wf(wf_vac)

    for basis in ccsdt_bases:
        base_name = mol_file.split(".")[0] + "_" + "CCSD(T)" + "_" + basis
        params = {"basis_set": basis,
                  "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                               "thresh": 14,
                                               "method": "ccsd(t)"}}}

        if mol._nelectrons < 100:
            params["overwrite_inputs"]["rem"]["mem_total"] = 760000
            params["overwrite_inputs"]["rem"]["cc_memory"] = 600000
            size = 1
        else:
            params["overwrite_inputs"]["rem"]["mem_total"] = 175000
            params["overwrite_inputs"]["rem"]["cc_memory"] = 140000
            size = 2

        name_vac = base_name + "_" + "vacuum"
        if name_vac not in finished:
            fw_vac = SinglePointFW(copy.deepcopy(mol),
                                   name=name_vac,
                                   qchem_input_params=copy.deepcopy(params),
                                   db_file=">>db_file<<")
            wf_vac = Workflow([fw_vac], name=name_vac)
            wf_vac = add_tags(wf_vac, {"class": "ssbt",
                                       "set": "20230702_mv_sp_wf",
                                       "size": size})
            print(mol._nelectrons, wf_vac)
            lp.add_wf(wf_vac)

    for basis in mp2_bases:
        base_name = mol_file.split(".")[0] + "_" + "MP2" + "_" + basis
        params = {"basis_set": basis,
                  "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                               "thresh": 14,
                                               "method": "mp2"}}}

        if mol._nelectrons < 100:
            params["overwrite_inputs"]["rem"]["mem_total"] = 760000
            params["overwrite_inputs"]["rem"]["cc_memory"] = 600000
            size = 1
        else:
            params["overwrite_inputs"]["rem"]["mem_total"] = 175000
            params["overwrite_inputs"]["rem"]["cc_memory"] = 140000
            size = 2

        name_vac = base_name + "_" + "vacuum"
        if name_vac not in finished:
            fw_vac = SinglePointFW(copy.deepcopy(mol),
                                   name=name_vac,
                                   qchem_input_params=copy.deepcopy(params),
                                   db_file=">>db_file<<")
            wf_vac = Workflow([fw_vac], name=name_vac)
            wf_vac = add_tags(wf_vac, {"class": "ssbt",
                                       "set": "20230702_mv_sp_wf",
                                       "size": size})
            print(mol._nelectrons, wf_vac)
            lp.add_wf(wf_vac)