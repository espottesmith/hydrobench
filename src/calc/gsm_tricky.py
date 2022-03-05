import copy

from monty.serialization import loadfn

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import (PESScanFW, FrequencyFlatteningTransitionStateFW,
                                          OptimizeFW)
from atomate.vasp.powerups import add_tags

from mpcat.automate.atomate.fireworks.core import SingleEndedGSMFW

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")

params = {"dft_rung": 4,
        "basis_set": "def2-svpd",
         "pcm_dielectric": 18.5,
         "overwrite_inputs": {"rem": {"method": "scan",
                                       "thresh": "14",
                                       "scf_algorithm": "diis"}}}

base_dir = "/Users/ewcss/data/ssbt/mol_files/"

isomers = {"acetalrct1_1": {"bonds_broken": [(22, 24)],
                            "bonds_formed": [(4, 24)]},
           "acetalpro1_1": {"bonds_broken": [(4, 24)],
                            "bonds_formed": [(22, 24)]},
           "acetalrct2_1": {"bonds_broken": [(1, 4)]},
           "acetalpro22_1": {"bonds_formed": [(1, 4)]},
           "acetalrct3_1": {"bonds_formed": [(1, 17)]},
           "acetalpro3_1": {"bonds_broken": [(1, 17)]},
           "acetalrct4_1": {"bonds_broken": [(17, 19)],
                            "bonds_formed": [(19, 20)]},
           "acetalpro4_1": {"bonds_broken": [(20, 21)],
                            "bonds_formed": [(17, 21)]},
           "acetalrct5_1": {"bonds_broken": [(20, 22)],
                            "bonds_formed": [(3, 22)]},
           "acetalpro5_1": {"bonds_broken": [(3, 22)],
                            "bonds_formed": [(20, 22)]},
           "acetalrct6_1": {"bonds_broken": [(1, 3)]},
           "acetalpro6_1": {"bonds_formed": [(1, 3)]},
           "ester1rct_-1": {"bonds_formed": [(6, 24)]},
           "ester1pro_-1": {"bonds_broken": [(6, 24)]},
           "ester2rct_-1": {"bonds_broken": [(6, 11)]},
           "ester2pro_-1": {"bonds_formed": [(6, 11)]}}

for filename, iso in isomers.items():
    mol = Molecule.from_file(base_dir + filename + ".xyz")
    mol.set_charge_and_spin(int(filename.split("_")[-1]))

    par = copy.deepcopy(params)

    name = filename + " GSM"

    fw = SingleEndedGSMFW(mol,
                          iso,
                            name=name,
                          max_cores=40,
                            input_params=par,
                            db_file=">>db_file<<"
                            )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20211002_gsm"})

    print(wf)
    lp.add_wf(wf)