import copy
import os

import numpy as np

from monty.serialization import dumpfn

from atomate.qchem.database import QChemCalcDb


def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)


db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

base_dir = "../data/dft/newmv/"

# Make naming consistent
replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

dft_data = [e for e in db.db["tasks"].find({"tags.set": "20220309_mv_sp_dft"}) if "SMD" not in e["task_label"]]
for d in dft_data:
    try:
        del d["custodian"]
    except:
        continue

dumpfn(dft_data, os.path.join(base_dir, "dft_sp.json"))