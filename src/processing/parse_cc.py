import os
import re

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

hf = re.compile(r"\s+SCF energy\s+=\s+([\-\.0-9]+)")
mp2 = re.compile(r"\s+MP2 energy\s+=\s+([\-\.0-9]+)")
ccsd_corr = re.compile(r"\s+CCSD correlation energy\s+=\s+([\-\.0-9]+)")
ccsd_total = re.compile(r"\s+CCSD total energy\s+=\s+([\-\.0-9]+)")
ccsdt_corr = re.compile(r"\s+CCSD\(T\) correlation energy\s+=\s+([\-\.0-9]+)")
ccsdt_total = re.compile(r"\s+CCSD\(T\) total energy\s+=\s+([\-\.0-9]+)")

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

all_data = dict()

for base_dir in ["20211016_sp_cc_dir",
                 "20211021_sp_cc_dir",
                 ]:
    for direct in os.listdir(base_dir):
        files = os.listdir(os.path.join(base_dir, direct))

        # Nothing ran
        if not any(["FW.json" in x for x in files]):
            continue
        # Error
        if "error.1.tar.gz" in files:
            continue

        if "FW.json" in files:
            meta = loadfn(os.path.join(base_dir, direct, "FW.json"))
        else:
            meta = loadfn(os.path.join(base_dir, direct, "FW.json.gz"))

        name = meta["name"]

        if "SMD" in name:
            continue

        all_data[name] = dict()

        if "mol.qout" in files:
            filename = "mol.qout"
        else:
            filename = "mol.qout.gz"

        with zopen(os.path.join(base_dir, direct, filename), mode="rt", encoding="ISO-8859-1") as f:
            text = f.read()

            thishf = hf.search(text)
            if thishf is not None:
                all_data[name]["hf"] = float(thishf.group(1))

            thismp2 = mp2.search(text)
            if thismp2 is not None:
                all_data[name]["mp2"] = float(thismp2.group(1))

            thisccsd_corr = ccsd_corr.search(text)
            if thisccsd_corr is not None:
                all_data[name]["ccsd_corr"] = float(thisccsd_corr.group(1))

            thisccsd_total = ccsd_total.search(text)
            if thisccsd_total is not None:
                all_data[name]["ccsd_total"] = float(thisccsd_total.group(1))

            thisccsdt_corr = ccsdt_corr.search(text)
            if thisccsdt_corr is not None:
                all_data[name]["ccsdt_corr"] = float(thisccsdt_corr.group(1))

            thisccsdt_total = ccsdt_total.search(text)
            if thisccsdt_total is not None:
                all_data[name]["ccsdt_total"] = float(thisccsdt_total.group(1))

dumpfn(all_data, "../data/reference/all_cc_data.json")