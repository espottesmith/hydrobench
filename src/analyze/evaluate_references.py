import copy
import statistics

from monty.serialization import loadfn, dumpfn

from pymatgen.core.structure import Molecule

data = loadfn("../data/dft/new/dft_compiled_data.json")

for_barriers = list()
rev_barriers = list()
dgs = list()
abs_dgs = list()

for rxn in data:
    try:
        for_barrier = data[rxn]["vacuum"]["reference"]["ts"] - data[rxn]["vacuum"]["reference"]["rct"]
        rev_barrier = data[rxn]["vacuum"]["reference"]["ts"] - data[rxn]["vacuum"]["reference"]["pro"]
        dg = data[rxn]["vacuum"]["reference"]["pro"] - data[rxn]["vacuum"]["reference"]["rct"]
        print(rxn, for_barrier, rev_barrier, dg)

        for_barriers.append(for_barrier)
        rev_barriers.append(rev_barrier)
        dgs.append(dg)
        abs_dgs.append(abs(dg))
    except TypeError:
        continue

barriers = for_barriers + rev_barriers

print("AVERAGES")
print("FORWARD BARRIERS", statistics.mean(for_barriers), statistics.stdev(for_barriers))
print("REVERSE BARRIERS", statistics.mean(rev_barriers), statistics.stdev(rev_barriers))
print("ALL BARRIERS", statistics.mean(barriers), statistics.stdev(barriers))
print("REACTION DG", statistics.mean(dgs), statistics.stdev(dgs))
print("|DG|", statistics.mean(abs_dgs), statistics.stdev(abs_dgs))