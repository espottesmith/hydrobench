import copy
import os

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import SinglePointFW

from atomate.common.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")
db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")


mp2_bases = ["def2-TZVPP", "def2-QZVPP"]
# ccsdt_bases = ["def2-TZVP"]

for calc in db.db["tasks"].find({"tags.set": {"$in": ["20220219_mp2_tsopt", "20220219_mp2_qirc"]}}):
    mol = Molecule.from_dict(calc["output"]["optimized_molecule"])

    name = calc["task_label"]

    if "ester" in name:
        continue

    for basis in mp2_bases:
        this_name = name.replace("def2-SVP", basis) + "SP"
        params = {"dft_rung": 4,
                  "basis_set": basis,
                  "overwrite_inputs": {"rem": {"method": "mp2",
                                               "thresh": "14",
                                               "scf_algorithm": "diis"}}}

        fw = SinglePointFW(copy.deepcopy(mol),
                           this_name,
                           qchem_input_params=params,
                           db_file=">>db_file<<")

        wf = Workflow([fw], name=this_name)
        wf = add_tags(wf, {"class": "ssbt", "set": "20220228_mp2_sp"})

        print(wf)
        lp.add_wf(wf)

    # for basis in ccsdt_bases:
    #     this_name = name.replace("def2-SVP", basis).replace("MP2", "CCSD(T)") + "SP"
    #     params = {"dft_rung": 4,
    #               "basis_set": basis,
    #               "overwrite_inputs": {"rem": {"method": "ccsd(t)",
    #                                            "thresh": "14",
    #                                            "scf_algorithm": "diis"}}}
    #
    #     fw = SinglePointFW(copy.deepcopy(mol),
    #                        this_name,
    #                        qchem_input_params=params,
    #                        db_file=">>db_file<<")
    #
    #     wf = Workflow([fw], name=this_name)
    #     wf = add_tags(wf, {"class": "ssbt", "set": "20220228_mp2_sp"})
    #
    #     print(wf)
    #     lp.add_wf(wf)