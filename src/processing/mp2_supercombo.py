from monty.serialization import loadfn, dumpfn


mp2 = loadfn("../data/reference/mp2_corrected_full.json")
cc = loadfn("../data/reference/grouped_cc_data.json")

new_data = dict()

for k, v in mp2.items():
    if k in cc:
        if "CCSD(T)_def2-TZVP" in cc[k]:
            energy = v + (cc[k]["CCSD(T)_def2-TZVP"]["ccsdt_total"] - cc[k]["CCSD(T)_def2-TZVP"]["mp2"])
            new_data[k] = energy

dumpfn(new_data, "../data/reference/mp2_super.json")