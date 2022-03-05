import copy

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb

from atomate.qchem.fireworks.core import (OptimizeFW)

from atomate.vasp.powerups import add_tags

lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")

params = {"dft_rung": 4,
          "basis_set": "def2-svpd",
          "overwrite_inputs": {"rem": {"method": "scan",
                                       "thresh": "14",
                                       "scf_algorithm": "diis"}}}


# base_dir = "/Users/ewcss/data/ssbt/mol_files/"
base_dir = "/Users/ewcss/data/ssbt/20220112_new_ts_guesses/"

# constraints = {"aceticanhydride_0": ["6 7 1.06", "4 7 1.42"],
#                "cl2co_0": ["1 4 2.386", "1 5 2.948", "5 7 0.996", "7 8 1.613"],
#                "cytosine1_0": ["10 14 1.53", "14 20 1.271", "17 20 1.151", "11 19 1.999"],
#                "cytosine2_0": ["10 11 1.922"],
#                "diazoniumsn2_1": ["4 12 2.230", "4 14 2.233"],
#                "furannoh2o1_1": ["16 18 1.371", "3 16 1.291"],
#                "furannoh2o2_1": ["1 3 2.030", "3 17 2.042", "17 19 0.983", "19 20 1.750"],
#                "furannoh2o3_1": ["4 16 1.211", "16 18 1.565"],
#                "furannoh2o4_1": ["3 17 2.210", "17 19 0.979", "19 20 1.780"],
#                "furannoh2o5_1": ["5 21 1.346", "20 21 1.315"],
#                "furannoh2o6_1": ["17 18 1.014", "18 19 1.551"],
#                "lactone_0": ["1 7 2.006", "1 10 2.024"]}

# constraints = {"aceticanhydride_0": ["6 7 1.06", "4 7 1.42", "1 2 1.820", "1 6 1.764"]}

# constraints = {"Cl2CO_TS1": ["2 6 1.331", "1 5 1.625", "5 6 1.200"],
#                "epoxide_ts1": ["1 3 1.940", "1 10 2.159"],
#                "epoxide_ts2": ["1 3 1.879", "1 10 1.982"]}

# constraints = {"enamine_ts_alex_0": ["4 14 1.725", "14 16 1.070", "13 16 1.355"]}

# constraints = {"lactone_basicmod1_-1": ["1 4 1.854"],
#                "lactone_basicmod2_-1": ["1 4 1.870"],
#                "lactone_basicmod3_-1": ["1 3 1.372", "1 8 1.832"],
#                "lactone_basicmod4_-1": ["1 9 1.764"],
#                "lactone_basicmod5_-1": ["1 3 1.386", "1 9 1.743"],
#                "lactone_basicmod6_-1": ["1 3 1.392", "1 13 1.709"],
#                "lactone_basicmod7_-1": ["3 7 1.876"]}

# constraints = {"imine1_0": ["2 8 1.354", "3 4 1.836", "4 8 1.758"],
#                "iminium2_1": ["3 4 1.503", "2 8 1.463", "4 8 1.155"],
#                "methylformate_-1": ["1 4 1.897", "1 9 1.381"]}

constraints = {"amide_aa_1": ["8 9 2.094", "9 18 1.681", "10 22 2.415"],
               "amide_aa_2": ["7 14 1.970", "7 23 1.959", "22 23 0.988", "18 22 2.228"],
               "amide_maa_1": ["2 13 1.803", "13 17 1.836", "3 16 2.141"],
               "amide_maa_2": ["2 9 1.937", "9 12 1.265", "12 16 1.248", "8 16 1.819"],
               "amide_mba_1": ["7 9 1.798", "9 24 1.846", "8 23 2.073"],
               "amide_mba_2": ["7 8 2.066", "8 24 1.355", "12 24 1.168", "12 18 1.767"],
               "carbonate_1": ["1 5 1.680", "1 7 2.173"],
               "borohydride_1": ["1 3 2.1317", "3 4 0.768", "4 5 2.821", "2 5 2.672"],
               "borohydride_2": ["1 5 2.425", "5 7 1.660"],
               "basic_epoxide_1": ["1 3 1.807", "1 8 2.102"],
               "basic_epoxide_2": ["1 3 1.881", "1 8 2.315"]}

# constraints = {"MeOAc_ts1_ac": ["3 8 1.721", "8 11 1.023", "9 11 1.546"],
#                "MeOAc_ts1_bc": ["1 7 2.471"]}

for filename, cons in constraints.items():
    mol = Molecule.from_file(base_dir + filename + ".xyz")
    # mol.set_charge_and_spin(int(filename.replace("ac", "1").replace("bc", "-1").split("_")[-1]))
    mol.set_charge_and_spin(-1)

    par = copy.deepcopy(params)
    parts = copy.deepcopy(params)
    par["overwrite_inputs"]["opt"] = {"CONSTRAINT": ["stre " + x for x in cons]}

    name = filename + " constrained: " + "; ".join(cons)
    namets = filename + " fftsopt"

    fw = OptimizeFW(molecule=copy.deepcopy(mol),
                    name=name,
                    qchem_input_params=par,
                    db_file=">>db_file<<"
                    )
    wf = Workflow([fw], name=name)
    wf = add_tags(wf, {"class": "ssbt", "set": "20220113_constrained_tsopt"})

    # fwts = FrequencyFlatteningTransitionStateFW(molecule=copy.deepcopy(mol),
    #                                             name=namets,
    #                                             qchem_input_params=parts,
    #                                             db_file=">>db_file<<"
    #                                             )
    # wfts = Workflow([fwts], name=namets)
    # wfts = add_tags(wfts, {"class": "ssbt", "set": "20211015_fftsopt"})

    print(wf)
    # print(wfts)
    lp.add_wf(wf)
    # lp.add_wf(wfts)