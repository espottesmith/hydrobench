import copy
import statistics

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule


all_sp = loadfn("../../data/dft/corrected/20230625_dft_sp.json")

reference_vac_pcm = loadfn("../../data/reference/corrected/20230625_mp2_super.json")

replacements = {"amide_1_1": "amide11",
                "amide_1_2": "amide12",
                "amide_2_1": "amide21",
                "amide_2_2": "amide22",
                "basic_epoxide_1": "basicepoxide1",
                "basic_epoxide_2": "basicepoxide2",
                }

compiled = dict()

rxns = set()
methods = set()


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

    if solv.lower() != "vacuum":
        continue

    rxns.add(rxn)
    methods.add(method)

    if rxn not in compiled:
        compiled[rxn] = {"vacuum": dict()}

    if method not in compiled[rxn][solv]:
        compiled[rxn][solv][method] = {"rct": None,
                                       "ts": None,
                                       "pro": None}

    compiled[rxn][solv][method][state] = energy

for rxn in rxns:
    for solv in ["vacuum"]:
        compiled[rxn][solv]["reference"] = {"rct": None,
                                            "ts": None,
                                            "pro": None}
        if solv == "SMD":
            continue
        else:
            compiled[rxn][solv]["reference"]["rct"] = reference_vac_pcm.get("{}_rct_{}".format(rxn, solv))
            compiled[rxn][solv]["reference"]["ts"] = reference_vac_pcm.get("{}_ts_{}".format(rxn, solv))
            compiled[rxn][solv]["reference"]["pro"] = reference_vac_pcm.get("{}_pro_{}".format(rxn, solv))

dumpfn(compiled, "../../data/dft/corrected/20230819_dg_compiled_data.json")

rxns = sorted(list(rxns))
methods = list(methods)

order = list()

references = dict()
for rxn in rxns:
    for solv in ["vacuum"]:
        name = "_".join([rxn, solv])
        order.append(name)

methods_errs = dict()
methods_abserrs = dict()
for method in methods:
    methods_errs[method] = list()
    methods_abserrs[method] = list()
    for rxn in rxns:
        for solv in ["vacuum"]:
            name = "_".join([rxn, solv])

            try:
                forward = compiled[rxn][solv][method]["pro"] - compiled[rxn][solv][method]["rct"]
            except (TypeError, KeyError) as e:
                forward = None

            try:
                forref = compiled[rxn][solv]["reference"]["pro"] - compiled[rxn][solv]["reference"]["rct"]
            except (TypeError, KeyError) as e:
                forref = None

            if forref is not None and forward is not None:
                forerr = (forward - forref) / forref
                methods_errs[method].append(forerr)
                methods_abserrs[method].append(abs(forerr))
            else:
                methods_errs[method].append(None)
                methods_abserrs[method].append(None)

# print(order)
# print(methods_abserrs["r2SCAN3c"])

for solv in ["vacuum"]:
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
    with open("../../data/results/corrected/dg_errs_rel_{}.csv".format(solv), "w") as file:
        file.write(",".join(header) + "\n")
        for line in alldata:
            file.write(",".join(line) + "\n")

for solv in ["vacuum"]:
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
    with open("../../data/results/corrected/dg_abserrs_rel_{}.csv".format(solv), "w") as file:
        file.write(",".join(header) + "\n")
        for line in alldata:
            file.write(",".join(line) + "\n")