import os
import re
import copy

from monty.io import zopen
from monty.serialization import loadfn, dumpfn

from .utils import (correct_hf, correct_wft,
                    alpha_tzvpp_qzvpp,
                    beta_tzvpp_qzvpp)


new_mp2 = loadfn("../data/reference/corrected/20230622_dumped_mp2.json")
new_cc = loadfn("../data/reference/corrected/20230622_dumped_cc.json")

groups = dict()

for n, data in all_data.items():

    name = copy.deepcopy(n)

    for a, b in replacements.items():
        name = name.replace(a, b)

    contents = name.split("_")
    groupname = "_".join([contents[0], contents[1], contents[-1]])
    try:
        basis = contents[4]
    except IndexError:
        print(n)
        continue

    if groupname not in groups:
        groups[groupname] = {basis: data}
    else:
        groups[groupname][basis] = data

dumpfn(groups, "grouped_mp2_data.json")

full_corrected_data = dict()

for groupname, groupdata in groups.items():
    try:
        hf = correct_hf(groupdata["def2-TZVPP"]["hf"],
                        groupdata["def2-QZVPP"]["hf"],
                        3,
                        4,
                        alpha=alpha_tzvpp_qzvpp)

        mp2 = correct_wft(groupdata["def2-TZVPP"]["mp2_corr"],
                           groupdata["def2-QZVPP"]["mp2_corr"],
                           3,
                           4,
                           beta=beta_tzvpp_qzvpp)

        total_e = hf + mp2
        full_corrected_data[groupname] = total_e
    except KeyError:
        print("PROBLEM CORRECTED FULL", groupname)
        continue

dumpfn(full_corrected_data, "mp2_corrected_full.json")