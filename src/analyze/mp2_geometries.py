import os
import copy

import numpy as np

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from atomate.qchem.database import QChemCalcDb


def perturb(mol, mode, scale=0.6):
    mol_copy = copy.deepcopy(mol)
    for ii in range(len(mol)):
        vec = np.array(mode[ii])
        mol_copy.translate_sites(indices=[ii], vector=vec * scale)
    return mol_copy


db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

base_dir = "/Users/ewcss/data/ssbt/for_sp"

mp2_data = {"imine": dict(),
            "carbonate": dict(),
            "amide": dict(),
            "diazonium": dict()}

baseline_mg = {
    "amide": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "amide_2_2_pro_-1.xyz")), OpenBabelNN())},
    "carbonate": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_rct_-1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_ts_-1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "carbonate_pro_-1.xyz")), OpenBabelNN())},
    "diazonium": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_rct_1.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_ts_1.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "diazonium_pro_1.xyz")), OpenBabelNN())},
    # "ester": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_rct_-1.xyz")), OpenBabelNN()),
    #           "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_ts_-1.xyz")), OpenBabelNN()),
    #           "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "ester_pro_-1.xyz")), OpenBabelNN())},
    "imine": {"rct": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_rct_0.xyz")), OpenBabelNN()),
              "ts": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_ts_0.xyz")), OpenBabelNN()),
              "pro": MoleculeGraph.with_local_env_strategy(Molecule.from_file(os.path.join(base_dir, "imine_pro_0.xyz")), OpenBabelNN())}
}

for ts_calc in db.db["tasks"].find({"tags.set": "20220219_mp2_tsopt"}):
    if "ester" in ts_calc["task_label"]:
        continue

    this_rname = None
    for rname in mp2_data:
        if rname in ts_calc["task_label"]:
            this_rname = rname
            break

    if this_rname is None:
        raise ValueError("No fitting reaction for", calc["task_label"])

    ts_mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(ts_calc["output"]["optimized_molecule"]), OpenBabelNN())
    mode = ts_calc["output"]["frequency_modes"][0]

    for_mol = perturb(ts_mg.molecule, mode, scale=0.6)
    rev_mol = perturb(ts_mg.molecule, mode, scale=-0.6)

    for_calc = None
    for_mg = None
    rev_calc = None
    rev_mg = None

    #TODO: finish this