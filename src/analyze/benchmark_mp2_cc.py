import copy
import os
import statistics

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

all_sp = loadfn("../data/dft/ssbt_dft_sp.json")

t_corrected = loadfn("../data/reference/cc_corrected_t.json")
small_corrected = loadfn("../data/reference/cc_corrected_small.json")
full_corrected = loadfn("../data/reference/cc_corrected_full.json")

mp2 = loadfn("../data/reference/mp2_super.json")

compiled = dict()

rxns = set()
solvs = ["vacuum", "IEF-PCM"]

types = {"small": small_corrected,
         "t": t_corrected,
         "full": full_corrected,
         "mp2": mp2}

for datapoint in all_sp:
    name = datapoint["task_label"]
    contents = name.split("_")
    rxn = contents[0]
    state = contents[1]
    method = contents[3]
    solv = contents[-1]
    energy = datapoint["output"]["final_energy"]

    rxns.add(rxn)

for rxn in rxns:
    compiled[rxn] = dict()
    for solv in solvs:
        compiled[rxn][solv] = dict()
        for method, data in types.items():

            compiled[rxn][solv][method] = {"rct": None,
                                           "ts": None,
                                           "pro": None}

            compiled[rxn][solv][method]["rct"] = data.get("{}_rct_{}".format(rxn, solv))
            compiled[rxn][solv][method]["ts"] = data.get("{}_ts_{}".format(rxn, solv))
            compiled[rxn][solv][method]["pro"] = data.get("{}_pro_{}".format(rxn, solv))

dumpfn(compiled, "../data/reference/wavefunction_compiled.json")

rxns = sorted(list(rxns))

order = list()

data = dict()
for rxn in rxns:
    for solv in solvs:
        forname = "_".join([rxn, solv, "forward"])
        revname = "_".join([rxn, solv, "reverse"])
        data[forname] = dict()
        data[revname] = dict()
        for method in types:
            try:
                forward = compiled[rxn][solv][method]["ts"] - compiled[rxn][solv][method]["rct"]
                data[forname][method] = forward * 27.2114
            except TypeError:
                data[forname][method] = None

            try:
                reverse = compiled[rxn][solv][method]["ts"] - compiled[rxn][solv][method]["pro"]
                data[revname][method] = reverse * 27.2114
            except TypeError:
                data[revname][method] = None

errors_t = dict()
errors_small = dict()
errors_mp2 = dict()

for nn, dd in data.items():
    if any([x is None for x in dd.values()]):
        continue

    errors_t[nn] = dd["t"] - dd["full"]
    errors_small[nn] = dd["small"] - dd["full"]
    errors_mp2[nn] = dd["mp2"] - dd["full"]

avg_error_t = statistics.mean(errors_t.values())
avg_error_small = statistics.mean(errors_small.values())
avg_error_mp2 = statistics.mean(errors_mp2.values())
mae_t = statistics.mean([abs(i) for i in errors_t.values()])
mae_small = statistics.mean([abs(i) for i in errors_small.values()])
mae_mp2 = statistics.mean([abs(i) for i in errors_mp2.values()])

for nn in errors_t:
    print("{}: CCSD TZVPP, (T) SVP/TZVP corrected: {}; CCSD(T) SVP/TZVP corrected: {}; MP2: {}".format(nn, errors_t[nn], errors_small[nn], errors_mp2[nn]))

print("Avg Error: CCSD TZVPP, (T) SVP/TZVP corrected: {}; CCSD(T) SVP/TZVP corrected: {}; MP2: {}".format(avg_error_t, avg_error_small, avg_error_mp2))
print("Avg Unsigned Error: CCSD TZVPP, (T) SVP/TZVP corrected: {}; CCSD(T) SVP/TZVP corrected: {}; MP2: {}".format(mae_t, mae_small, mae_mp2))
