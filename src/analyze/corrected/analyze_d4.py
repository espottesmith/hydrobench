import copy
import statistics

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule

from atomate.qchem.database import QChemCalcDb


db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")

reference_vac_pcm = loadfn("../../data/reference/corrected/20230625_mp2_super.json")

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                "realts": "ts"
                }

compiled = dict()

rxns = set()
methods = set()
solvs = ["vacuum", "IEF-PCM"]

all_sp = db.db["tasks"].find({"tags.set": {"$in": ["20230403_sp_d4"]}})

for datapoint in all_sp:
    name = copy.deepcopy(datapoint["task_label"])
    for a, b in replacements.items():
        name = name.replace(a, b)

    if name.startswith("_"):
        name = "epoxide2_pro_1" + name

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

#override old aceticanhydride calcs
for datapoint in db.db["tasks"].find({"tags.set": {"$in": ["20230702_sp_d4"]}}):
    name = copy.deepcopy(datapoint["task_label"])
    for a, b in replacements.items():
        name = name.replace(a, b)

    if name.startswith("_"):
        name = "epoxide2_pro_1" + name

    print("REPLACEMENT", name)

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

dumpfn(compiled, "../../data/dft/corrected/d4_dft_compiled_data.json")

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
                forerr = (forward - forref) * 27.2114
                methods_errs[method].append(forerr)
                methods_abserrs[method].append(abs(forerr))
            else:
                methods_errs[method].append(None)
                methods_abserrs[method].append(None)

            if revref is not None and reverse is not None:
                reverr = (reverse - revref) * 27.2114
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
            print("STATS ERROR", method, solv)
            continue
        alldata.append(line)
    alldata = sorted(alldata, key=lambda x: float(x[-1]))
    with open("../../data/results/corrected/d4_errs_{}.csv".format(solv), "w") as file:
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
    with open("../../data/results/corrected/d4_abserrs_{}.csv".format(solv), "w") as file:
        file.write(",".join(header) + "\n")
        for line in alldata:
            file.write(",".join(line) + "\n")