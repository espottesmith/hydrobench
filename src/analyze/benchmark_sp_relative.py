import copy
import os
import statistics

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

all_sp = loadfn("/Users/ewcss/data/ssbt/20220204_sp_benchmark/dft_sp.json")

# reference_vac_pcm = loadfn("/Users/ewcss/data/ssbt/20220205_reparse/cc_corrected_full.json")
# reference_vac_pcm = loadfn("/Users/ewcss/data/ssbt/20220205_reparse/cc_corrected_t.json")
reference_vac_pcm = loadfn("/Users/ewcss/data/ssbt/20220211_mp2/mp2_super.json")

r2scan3c_tasks = loadfn("/Users/ewcss/data/ssbt/20220211_r2scan3c/r2scan3c_energies_defgrid3.json")

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

for calc in r2scan3c_tasks:
    orig_name = copy.deepcopy(calc["task_label"])
    name = calc["task_label"]
    for a, b in replacements.items():
        name = name.replace(a, b)

    calc["task_label"] = name + "_r2SCAN3c_vacuum"
    calc["output"]["final_energy"] /= 27.2114

    # print(orig_name, calc["task_label"])
    if calc["task_label"].startswith("_"):
        continue

    all_sp.append(calc)

compiled = dict()

rxns = set()
methods = set()
solvs = ["vacuum", "IEF-PCM"]

for datapoint in all_sp:
    name = datapoint["task_label"]
    contents = name.split("_")
    rxn = contents[0]
    state = contents[1]
    method = contents[3]
    solv = contents[-1]
    energy = datapoint["output"]["final_energy"]

    rxns.add(rxn)
    methods.add(method)

    if rxn not in compiled:
        compiled[rxn] = {"vacuum": dict(),
                         "IEF-PCM": dict()}

    if method not in compiled[rxn][solv]:
        compiled[rxn][solv][method] = {"rct": None,
                                       "ts": None,
                                       "pro": None}

    compiled[rxn][solv][method][state] = energy

for rxn in rxns:
    for solv in solvs:
        compiled[rxn][solv]["reference"] = {"rct": None,
                                            "ts": None,
                                            "pro": None}
        if solv == "SMD":
            continue
        else:
            compiled[rxn][solv]["reference"]["rct"] = reference_vac_pcm.get("{}_rct_{}".format(rxn, solv))
            compiled[rxn][solv]["reference"]["ts"] = reference_vac_pcm.get("{}_ts_{}".format(rxn, solv))
            compiled[rxn][solv]["reference"]["pro"] = reference_vac_pcm.get("{}_pro_{}".format(rxn, solv))

dumpfn(compiled, "/Users/ewcss/data/ssbt/20220211_benchmark/dft_compiled_data.json")

rxns = sorted(list(rxns))
methods = list(methods)

order = list()

references = dict()
for rxn in rxns:
    for solv in solvs:
        forname = "_".join([rxn, solv, "forward"])
        revname = "_".join([rxn, solv, "reverse"])
        order.append(forname)
        order.append(revname)

methods_errs = dict()
methods_abserrs = dict()
for method in methods:
    methods_errs[method] = list()
    methods_abserrs[method] = list()
    for rxn in rxns:
        for solv in solvs:
            forname = "_".join([rxn, solv, "forward"])
            revname = "_".join([rxn, solv, "reverse"])

            try:
                forward = compiled[rxn][solv][method]["ts"] - compiled[rxn][solv][method]["rct"]
            except (TypeError, KeyError) as e:
                forward = None

            try:
                reverse = compiled[rxn][solv][method]["ts"] - compiled[rxn][solv][method]["pro"]
            except (TypeError, KeyError) as e:
                reverse = None

            try:
                forref = compiled[rxn][solv]["reference"]["ts"] - compiled[rxn][solv]["reference"]["rct"]
            except (TypeError, KeyError) as e:
                forref = None

            try:
                revref = compiled[rxn][solv]["reference"]["ts"] - compiled[rxn][solv]["reference"]["pro"]
            except (TypeError, KeyError) as e:
                revref = None

            if forref is not None and forward is not None:
                forerr = (forward - forref) / forref
                methods_errs[method].append(forerr)
                methods_abserrs[method].append(abs(forerr))
            else:
                methods_errs[method].append(None)
                methods_abserrs[method].append(None)

            if revref is not None and reverse is not None:
                reverr = (reverse - revref) / revref
                methods_errs[method].append(reverr)
                methods_abserrs[method].append(abs(reverr))
            else:
                methods_errs[method].append(None)
                methods_abserrs[method].append(None)

# print(order)
# print(methods_abserrs["r2SCAN3c"])

for solv in solvs:
    appro_indices = [i for i, e in enumerate(order) if solv in e]
    header = ["method"] + [e for i, e in enumerate(order) if i in appro_indices] + ["Average"]
    alldata = list()
    for method, data in methods_errs.items():
        this_data = [d for i, d in enumerate(data) if i in appro_indices]
        try:
            line = [method] + [str(d) for d in this_data] + [str(statistics.mean([d for d in this_data if d is not None]))]
        except statistics.StatisticsError:
            continue
        alldata.append(line)
    alldata = sorted(alldata, key=lambda x: float(x[-1]))
    with open("/Users/ewcss/data/ssbt/20220211_benchmark/errs_rel_{}.csv".format(solv), "w") as file:
        file.write(",".join(header) + "\n")
        for line in alldata:
            file.write(",".join(line) + "\n")

for solv in solvs:
    appro_indices = [i for i, e in enumerate(order) if solv in e]
    header = ["method"] + [e for i, e in enumerate(order) if i in appro_indices] + ["Average"]
    alldata = list()
    for method, data in methods_abserrs.items():
        this_data = [d for i, d in enumerate(data) if i in appro_indices]
        try:
            line = [method] + [str(d) for d in this_data] + [str(statistics.mean([d for d in this_data if d is not None]))]
        except statistics.StatisticsError:
            print("STATS ERROR", method, solv)
            continue
        alldata.append(line)
    alldata = sorted(alldata, key=lambda x: float(x[-1]))
    with open("/Users/ewcss/data/ssbt/20220211_benchmark/abserrs_rel_{}.csv".format(solv), "w") as file:
        file.write(",".join(header) + "\n")
        for line in alldata:
            file.write(",".join(line) + "\n")