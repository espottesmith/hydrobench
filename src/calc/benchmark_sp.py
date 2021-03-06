import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

# Geometry optimization: SCAN/def2-SVPD

# General parameters:
# - THRESH 14 (except CC?)

# For CC:
# CCSD: def2-TZVPP, def2-QZVPP
# CCSD(T): def2-SVP, def2-TZVP

# GGA (10):
# - PBE
# - PBE-D3
# - BLYP
# - BLYP-D3
# - B97-D
# - B97-D3(0)
# - mPW91
# - mPW91-D3
# - VV10
# - rVV10
#
# Meta-GGA (7):
# - M06-L
# - M06-L-D3
# - SCAN
# - SCAN-D3
# - TPSS
# - TPSS-D3
# - MN12-L
# - B97M-rV
#
# Hybrid GGA (11):
# - PBE0 (global)
# - PBE0-D3 (global)
# - B3LYP (global)
# - B3LYP-D3 (global)
# - CAM-B3LYP (range-separated)
# - CAM-B3LYP-D3 (range-separated)
# - mPW1PW91 (global)
# - mPW1PW91-D3 (global)
# - wB97X (range-separated)
# - wB97X-D (range-separated)
# - wB97X-D3 (range-separated)
# - wB97X-V (range-separated)
#
# Hybrid Meta-GGA (13):
# - M06-2X (global)
# - M06-2X-D3 (global)
# - M06-HF (global)
# - M08-SO (global)
# - M11 (range-separated)
# - MN15 (global)
# - BMK (global)
# - BMK-D3 (global)
# - TPSSh (global)
# - SCAN0 (global)
# - SCAN0-D3 (global)
# - MPWB1K (global)
# - MPWB1K-D3 (global)
# - wB97M-V (range-separated)


dft_methods = ["PBE", ("PBE", "D3(BJ)"), "BLYP", ("BLYP", "D3(BJ)"), "B97-D", "B97-D3", "mPW91", ("mPW91", "D3(BJ)"), "VV10", "rVV10",
               "M06-L", ("M06-L", "D3(0)"), "SCAN", ("SCAN", "D3(BJ)"), "TPSS", ("TPSS", "D3(BJ)"), "MN12-L", ("MN12-L", "D3(BJ)"), "B97M-rV",
               "PBE0", ("PBE0", "D3(BJ)"), "B3LYP", ("B3LYP", "D3(BJ)"), "CAM-B3LYP", ("CAM-B3LYP", "D3(0)"), "mPW1PW91", ("mPW1PW91", "D3(BJ)"), "wB97X", "wB97XD", "wB97XD3", "wB97XV",
               "M06-2X", ("M06-2X", "D3(0)"), "M06-HF", "M08-SO", "M11", "MN15", "BMK", ("BMK", "D3(BJ)"), "TPSSh", ("TPSSh", "D3(BJ)"), "SCAN0", "mPWB1K", ("mPWB1K", "D3(BJ)"), "wB97M-V"]

dft_basis = "def2-TZVPPD"

ccsd_bases = ["def2-TZVPP", "def2-QZVPP"]
ccsdt_bases = ["def2-SVP", "def2-TZVP"]

mp2_bases = ["def2-TZVPP", "def2-QZVPP"]

base_dir = "/Users/ewcss/data/ssbt/for_sp"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20211016_sp_dft", "20211016_sp_cc"]}})]

for mol_file in os.listdir(base_dir):
    if "enamine" in mol_file:
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
                                   "set": "20211016_sp_dft"})
            print(wf_vac)
            lp.add_wf(wf_vac)

        name_pcm = base_name + "_" + "IEF-PCM"
        if name_pcm not in finished:
            pcm_params = copy.deepcopy(params)
            pcm_params["pcm_dielectric"] = 78.39
            pcm_params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
            fw_pcm = SinglePointFW(copy.deepcopy(mol),
                                   name=name_pcm,
                                   qchem_input_params=pcm_params,
                                   db_file=">>db_file<<")
            wf_pcm = Workflow([fw_pcm], name=name_pcm)
            wf_pcm = add_tags(wf_pcm, {"class": "ssbt",
                                       "set": "20211016_sp_dft"})
            print(wf_pcm)
            lp.add_wf(wf_pcm)

    for basis in ccsd_bases:
        base_name = mol_file.split(".")[0] + "_" + "CCSD" + "_" + basis
        params = {"basis_set": basis,
                  "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                               "thresh": 14,
                                               "method": "ccsd"}}}

        if mol._nelectrons < 100:
            params["overwrite_inputs"]["rem"]["mem_total"] = 750000
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
                                   "set": "20211016_sp_cc",
                                       "size": size})
            print(mol._nelectrons, wf_vac)
            lp.add_wf(wf_vac)

        name_pcm = base_name + "_" + "IEF-PCM"
        if name_pcm not in finished:
            pcm_params = copy.deepcopy(params)
            pcm_params["pcm_dielectric"] = 78.39
            pcm_params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
            fw_pcm = SinglePointFW(copy.deepcopy(mol),
                                   name=name_pcm,
                                   qchem_input_params=pcm_params,
                                   db_file=">>db_file<<")
            wf_pcm = Workflow([fw_pcm], name=name_pcm)
            wf_pcm = add_tags(wf_pcm, {"class": "ssbt",
                                       "set": "20211016_sp_cc",
                                       "size": size})
            print(mol._nelectrons, wf_pcm)
            lp.add_wf(wf_pcm)

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
                                   "set": "20211016_sp_cc",
                                       "size": size})
            print(mol._nelectrons, wf_vac)
            lp.add_wf(wf_vac)

        name_pcm = base_name + "_" + "IEF-PCM"
        if name_pcm not in finished:
            pcm_params = copy.deepcopy(params)
            pcm_params["pcm_dielectric"] = 78.39
            pcm_params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
            fw_pcm = SinglePointFW(copy.deepcopy(mol),
                                   name=name_pcm,
                                   qchem_input_params=pcm_params,
                                   db_file=">>db_file<<")
            wf_pcm = Workflow([fw_pcm], name=name_pcm)
            wf_pcm = add_tags(wf_pcm, {"class": "ssbt",
                                       "set": "20211016_sp_cc",
                                       "size": size})
            print(mol._nelectrons, wf_pcm)
            lp.add_wf(wf_pcm)

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
                                       "set": "20220209_mp2",
                                       "size": size})
            print(mol._nelectrons, wf_vac)
            lp.add_wf(wf_vac)

        name_pcm = base_name + "_" + "IEF-PCM"
        if name_pcm not in finished:
            pcm_params = copy.deepcopy(params)
            pcm_params["pcm_dielectric"] = 78.39
            pcm_params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
            fw_pcm = SinglePointFW(copy.deepcopy(mol),
                                   name=name_pcm,
                                   qchem_input_params=pcm_params,
                                   db_file=">>db_file<<")
            wf_pcm = Workflow([fw_pcm], name=name_pcm)
            wf_pcm = add_tags(wf_pcm, {"class": "ssbt",
                                       "set": "20220209_mp2",
                                       "size": size})
            print(mol._nelectrons, wf_pcm)
            lp.add_wf(wf_pcm)