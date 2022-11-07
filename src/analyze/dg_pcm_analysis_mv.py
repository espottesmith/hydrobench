import csv
import os

import difflib

base_dir = "../data/results/dg_mv"

vac_mae = dict()
vac_rel = dict()
pcm_mae = dict()
pcm_rel = dict()

with open(os.path.join(base_dir, "abserrs_IEF-PCM.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])
        pcm_mae[funct] = avg

with open(os.path.join(base_dir, "abserrs_rel_IEF-PCM.csv")) as file:
    reader = csv.reader(file)
    for i, row in enumerate(reader):
        if i == 0:
            continue
        elif row[0].lower() == "average" or "3c" in row[0].lower():
            continue
        funct = row[0]
        avg = float(row[-1])
        pcm_rel[funct] = avg

lcd = [e for e in pcm_rel if all([e in x for x in [pcm_mae]])]

for s in [pcm_mae, pcm_rel]:
    try:
        del s["Average"]
        del s["average"]
    except:
        continue

sorted_pcm_mae = sorted(pcm_mae.items(), key=lambda x: x[1])
sorted_pcm_rel = sorted(pcm_rel.items(), key=lambda x: x[1])

rankings = list()
for l in [sorted_pcm_mae, sorted_pcm_rel]:
    rankings.append([x[0] for x in l])


for r in rankings:
    print(r)

compiled = dict()

for dset in [sorted_pcm_mae, sorted_pcm_rel]:
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
        total.append(str(round(values[2 * i], 3)))
        total.append(str(values[2 * i + 1] + 1))
    final.append(total)

for t in final:
    print(" & ".join(t) + "\\\\")
    print("\\hline")