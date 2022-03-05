import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW, OptimizeFW

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

dft_methods = ["PBE", ("PBE", "D3(BJ)"), "BLYP", ("BLYP", "D3(BJ)"), "B97-D", "B97-D3", "mPW91", ("mPW91", "D3(BJ)"), "VV10", "rVV10",
               "M06-L", ("M06-L", "D3(0)"), "SCAN", ("SCAN", "D3(BJ)"), "TPSS", ("TPSS", "D3(BJ)"), "MN12-L", ("MN12-L", "D3(BJ)"), "B97M-rV",
               "PBE0", ("PBE0", "D3(BJ)"), "B3LYP", ("B3LYP", "D3(BJ)"), "CAM-B3LYP", ("CAM-B3LYP", "D3(0)"), "mPW1PW91", ("mPW1PW91", "D3(BJ)"), "wB97X", "wB97XD", "wB97XD3", "wB97XV",
               "M06-2X", ("M06-2X", "D3(0)"), "M06-HF", "M08-SO", "M11", "MN15", "BMK", ("BMK", "D3(BJ)"), "TPSSh", ("TPSSh", "D3(BJ)"), "SCAN0", "mPWB1K", ("mPWB1K", "D3(BJ)"), "wB97M-V"]

# Eventually reopt in def2-TZVPPD
dft_basis = "def2-SVPD"

constraints = {"diazonium_ts_1": ["4 12 2.255", "4 14 2.438"],
               "ester_ts_-1": ["1 7 2.628", "4 7 2.110", "3 4 1.110"],
               "imine_ts_0": ["2 8 1.062", "3 4 2.273", "4 8 1.705"],
               "carbonate_ts_-1": ["1 5 1.738", "1 7 2.110"],
               "amide_2_2_ts_-1": ["2 9 1.834", "4 8 1.005", "8 16 1.712", "9 12 1.208", "12 16 1.312"]}

base_dir = "/Users/ewcss/data/ssbt/for_sp"

finished = [e["task_label"] for e in db.db["tasks"].find({"tags.set": {"$in": ["20220130_opt_benchmark_constopts"]}})]

for mol_name, const in constraints.items():
    for mol_file in os.listdir(base_dir):
        if mol_name in mol_file:
            mol = Molecule.from_file(os.path.join(base_dir, mol_file))
            mol.set_charge_and_spin(int(mol_file.split(".")[0].split("_")[-1]))

            for method in dft_methods:
                params = {"basis_set": dft_basis,
                          "overwrite_inputs": {"rem": {"scf_algorithm": "diis",
                                                       "thresh": 14},
                                               "opt": {"CONSTRAINT": ["stre " + x for x in const]}}}
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
                    fw_vac = OptimizeFW(copy.deepcopy(mol),
                                           name=name_vac + "_constopt",
                                           qchem_input_params=copy.deepcopy(params),
                                           db_file=">>db_file<<")
                    wf_vac = Workflow([fw_vac], name=name_vac + "_constopt")
                    wf_vac = add_tags(wf_vac, {"class": "ssbt",
                                           "set": "20220130_opt_benchmark_constopts"})
                    print(wf_vac)
                    lp.add_wf(wf_vac)

                # name_pcm = base_name + "_" + "IEF-PCM"
                # if name_pcm not in finished:
                #     pcm_params = copy.deepcopy(params)
                #     pcm_params["pcm_dielectric"] = 78.39
                #     pcm_params["overwrite_inputs"]["pcm"] = {"theory": "IEFPCM"}
                #     fw_pcm = OptimizeFW(copy.deepcopy(mol),
                #                            name=name_pcm + "_constopt",
                #                            qchem_input_params=pcm_params,
                #                            db_file=">>db_file<<")
                #     wf_pcm = Workflow([fw_pcm], name=name_pcm + "_constopt")
                #     wf_pcm = add_tags(wf_pcm, {"class": "ssbt",
                #                                "set": "20220130_opt_benchmark_constopts"})
                #     print(wf_pcm)
                #     lp.add_wf(wf_pcm)