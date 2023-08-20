import copy
import os

import numpy as np

from monty.serialization import dumpfn

from atomate.qchem.database import QChemCalcDb


db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

base_dir = "../data/dft/corrected/"

# Make naming consistent
replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

og_data = [e for e in db.db["tasks"].find(
    {"tags.set": "20211016_sp_dft"}, 
    {"input": 1, "tags": 1, "task_label": 1, "output": 1, "formula_alphabetical": 1, "task_id": 1, "orig": 1}
) if "vacuum" in e["task_label"]]
new_data = [e for e in db.db["tasks"].find(
    {"tags.set": "20230617_sp_correct_dft"}, 
    {"input": 1, "tags": 1, "task_label": 1, "output": 1, "formula_alphabetical": 1, "task_id": 1, "orig": 1, "walltime": 1, "cputime": 1}
) if "vacuum" in e["task_label"]]

new_names = [x["task_label"].replace("realts", "ts") for x in new_data]

for ii, calc in enumerate(og_data):
    for jj, name in enumerate(new_names):
        if calc["task_label"] == name:
            print("REPLACING", ii, jj)
            og_data[ii] = new_data[jj]

dumpfn(og_data, os.path.join(base_dir, "20230625_dft_sp.json"))