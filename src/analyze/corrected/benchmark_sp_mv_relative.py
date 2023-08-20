import statistics

from monty.serialization import loadfn, dumpfn

all_sp = loadfn("../../data/dft/corrected/20230818_mv_dft_sp.json")
reference_vac_pcm = loadfn("../../data/reference/corrected/20230819_mvopt_svpd_reference_data.json")

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
    name = datapoint["task_label"]

    for a, b in replacements.items():
        name = name.replace(a, b)

    contents = name.split("_")
    rxn = contents[0]
    state = contents[1]
    method = contents[3]
    solv = contents[-1]
    energy = datapoint["output"]["final_energy"]

    rxns.add(rxn)
    methods.add(method)

    if solv != "vacuum":
        continue

    if rxn not in compiled:
        compiled[rxn] = dict()

    if method not in compiled[rxn]:
        compiled[rxn][method] = {"rct": None,
                                 "ts": None,
                                 "pro": None}

    compiled[rxn][method][state] = energy

# print("METHODS", methods)

for rxn in rxns:
    compiled[rxn]["reference"] = {"rct": None,
                                  "ts": None,
                                  "pro": None}

    compiled[rxn]["reference"]["rct"] = reference_vac_pcm[rxn]["rct"]["total_corrected"]
    compiled[rxn]["reference"]["ts"] = reference_vac_pcm[rxn]["ts"]["total_corrected"]
    compiled[rxn]["reference"]["pro"] = reference_vac_pcm[rxn]["pro"]["total_corrected"]

dumpfn(compiled, "../../data/dft/corrected/20230818_mv_dft_compiled_data.json")

rxns = sorted(list(rxns))
methods = list(methods)

order = list()

references = dict()
for rxn in rxns:
    forname = "_".join([rxn, "forward"])
    revname = "_".join([rxn, "reverse"])
    order.append(forname)
    order.append(revname)


# TODO: you are here
methods_errs = dict()
methods_abserrs = dict()
for method in methods:
    methods_errs[method] = list()
    methods_abserrs[method] = list()
    for rxn in rxns:
        forname = "_".join([rxn, "forward"])
        revname = "_".join([rxn, "reverse"])

        try:
            forward = compiled[rxn][method]["ts"] - compiled[rxn][method]["rct"]
        except (TypeError, KeyError) as e:
            problems = list()
            if compiled[rxn][method]["ts"] is None:
                problems.append("ts")
            if compiled[rxn][method]["rct"] is None:
                problems.append("rct")
            print("PROBLEM FORWARD", rxn, method, problems)
            forward = None

        try:
            reverse = compiled[rxn][method]["ts"] - compiled[rxn][method]["pro"]
        except (TypeError, KeyError) as e:
            problems = list()
            if compiled[rxn][method]["ts"] is None:
                problems.append("ts")
            if compiled[rxn][method]["pro"] is None:
                problems.append("pro")
            print("PROBLEM REVERSE", rxn, method, problems)
            reverse = None

        try:
            forref = compiled[rxn]["reference"]["ts"] - compiled[rxn]["reference"]["rct"]
        except (TypeError, KeyError) as e:
            print("PROBLEM FORWARD REFERENCE", rxn, method)
            forref = None

        try:
            revref = compiled[rxn]["reference"]["ts"] - compiled[rxn]["reference"]["pro"]
        except (TypeError, KeyError) as e:
            print("PROBLEM REVERSE REFERENCE", rxn, method)
            revref = None

        if forref is not None and forward is not None:
            forerr = (forward - forref) / forref
            methods_errs[method].append(forerr)
            methods_abserrs[method].append(abs(forerr))
        else:
            print("ADDING NONE FORWARD", forname)
            methods_errs[method].append(None)
            methods_abserrs[method].append(None)

        if revref is not None and reverse is not None:
            reverr = (reverse - revref) / revref
            methods_errs[method].append(reverr)
            methods_abserrs[method].append(abs(reverr))
        else:
            print("ADDING NONE REVERSE", revname)
            methods_errs[method].append(None)
            methods_abserrs[method].append(None)


header = ["method"] + [e for i, e in enumerate(order)] + ["Average"]
alldata = list()
for method, data in methods_errs.items():
    this_data = [d for i, d in enumerate(data)]
    try:
        line = [method] + [str(d) for d in this_data] + [str(statistics.mean([d for d in this_data if d is not None]))]
    except statistics.StatisticsError:
        continue
    alldata.append(line)
alldata = sorted(alldata, key=lambda x: float(x[-1]))
with open("../../data/results/corrected/mv_errs_relative_{}.csv".format(solv), "w") as file:
    file.write(",".join(header) + "\n")
    for line in alldata:
        file.write(",".join(line) + "\n")

header = ["method"] + [e for i, e in enumerate(order)] + ["Average"]
alldata = list()
for method, data in methods_abserrs.items():
    this_data = [d for i, d in enumerate(data)]
    try:
        line = [method] + [str(d) for d in this_data] + [str(statistics.mean([d for d in this_data if d is not None]))]
    except statistics.StatisticsError:
        print("STATS ERROR", method, solv)
        continue
    alldata.append(line)
alldata = sorted(alldata, key=lambda x: float(x[-1]))
with open("../../data/results/corrected/mv_abserrs_relative_{}.csv".format(solv), "w") as file:
    file.write(",".join(header) + "\n")
    for line in alldata:
        file.write(",".join(line) + "\n")
