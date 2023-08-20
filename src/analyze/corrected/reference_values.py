from monty.serialization import loadfn

data = loadfn("../../data/dft/corrected/20230625_dft_compiled_data.json")

sorted_rxns = sorted(list(data.keys()))

print("vacuum")
for r in sorted_rxns:
    d = data[r]
    if not any([x is None for x in d["vacuum"].get("reference", {"x": None}).values()]):
        for_barrier = (d["vacuum"]["reference"]["ts"] - d["vacuum"]["reference"]["rct"]) * 27.2114
        rev_barrier = (d["vacuum"]["reference"]["ts"] - d["vacuum"]["reference"]["pro"]) * 27.2114
        dg = (d["vacuum"]["reference"]["pro"] - d["vacuum"]["reference"]["rct"]) * 27.2114
        print(f"{r} & forward & {for_barrier:.3f} & {dg:.3f}\\\\\n\\hline")
        print(f"{r} & reverse & {rev_barrier:.3f} & {-1 * dg:.3f}\\\\\n\\hline")

print("\n\n\n")