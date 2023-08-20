import os

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

dft_data = [e for e in db.db["tasks"].find({"tags.set": {"$in": ["20220309_mv_sp_dft", "20230702_mv_sp_dft"]}}) if "SMD" not in e["task_label"]]
for d in dft_data:
    try:
        del d["custodian"]
    except:
        continue

dumpfn(dft_data, os.path.join(base_dir, "20230818_mv_dft_sp.json"))
