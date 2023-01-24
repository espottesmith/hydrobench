import csv
import os

base_dir = "../data/results/new"

vac_mae = dict()
vac_rel = dict()
pcm_mae = dict()
pcm_rel = dict()

with open(os.path.join(base_dir, "abserrs_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])
        vac_mae[funct] = avg

with open(os.path.join(base_dir, "abserrs_rel_vacuum.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])
        vac_rel[funct] = avg

lcd = [e for e in vac_mae if all([e in x for x in [vac_rel]])]

for s in [vac_mae, vac_rel]:
    try:
        del s["Average"]
        del s["average"]
    except:
        continue

sorted_vac_mae = sorted(vac_mae.items(), key=lambda x: x[1])
sorted_vac_rel = sorted(vac_rel.items(), key=lambda x: x[1])

compiled = dict()

for dset in [sorted_vac_mae, sorted_vac_rel]:
    for index, (name, val) in enumerate(dset):
        if name not in lcd:
            continue

        if name not in compiled:
            compiled[name] = list()

        compiled[name].append(val)
        compiled[name].append(index)

sorted_compiled = sorted(compiled.items(),
                         key=lambda x: sum([x[1][1], x[1][3]]) / 2)

final = list()
for name, values in sorted_compiled:
    total = [name]
    for i in range(2):
        total.append("{:.3f}".format(values[2 * i]))
        total.append(str(values[2 * i + 1] + 1))
    final.append(total)

for t in final:
    print(" & ".join(t) + "\\\\")
    print("\\hline")