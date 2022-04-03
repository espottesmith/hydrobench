import copy
import json
import os

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from monty.serialization import dumpfn, loadfn
from pandas import DataFrame as df
from pymongo import MongoClient

from fireworks.core.launchpad import LaunchPad
from pymatgen.analysis.molecule_matcher import (
    MoleculeMatcher,
    KabschMatcher,
    BruteForceOrderMatcher,
)
from pymatgen.core.structure import Molecule

# from pymatgen.core.units import Ha_to_eV

pio.renderers.default = "browser"
Ha_to_eV = 27.211
kcal_to_Ha = 1 / Ha_to_eV
data_dir = "../../../20220301_opt_benchmark"


def plot_heatmap(data, title="", zmin=None, zmax=None):
    fig = go.Figure()
    dataframe = df.from_dict(data)
    if not zmax:
        zmax = max(dataframe["average"].to_list())
    if not zmin:
        zmin = min(dataframe["average"].to_list())
    # max_val = dataframe.to_numpy().max()
    dataframe = dataframe.sort_values(by="average")
    fig.add_trace(
        go.Heatmap(
            z=dataframe,
            y=dataframe.index,
            x=dataframe.columns,
            colorscale="viridis",
            reversescale=True,
            zmin=zmin,
            zmax=zmax,
        )
    )

    fig.update_yaxes(
        autorange="reversed", visible=True, showgrid=False, nticks=len(dataframe.index)
    )
    fig.update_xaxes(side="top", visible=True, showgrid=False)
    fig.update_layout(width=800, height=800, plot_bgcolor="#808080", title=title)
    # fig.write_image("mae_heatmap_kcal.pdf")
    fig.show()


def get_elec_mae(reference, target):
    """

    :param reference:
    :param target:
    :return: MAE: Mean Absolute Error for electronic energies of barrier
    heights as a dict of {"molecule":{"functional": mae, . . . }, . . .,
    "average": {"functional": average_mae, . . .}}
    """
    mae = {}
    ref_barriers = {}
    target_barriers = {}
    for moltype, methods in target.items():

        mae[moltype + "_fwd"] = {}
        mae[moltype + "_rev"] = {}

        target_barriers[moltype + "_fwd"] = {}
        target_barriers[moltype + "_rev"] = {}

        ref_ts_e = reference[moltype]["ts"]["total_corrected"]
        ref_rct_e = reference[moltype]["rct"]["total_corrected"]
        ref_pro_e = reference[moltype]["pro"]["total_corrected"]
        ref_fwd_e = ref_ts_e - ref_rct_e
        ref_rev_e = ref_ts_e - ref_pro_e

        ref_barriers[moltype + "_fwd"] = ref_fwd_e * Ha_to_eV
        ref_barriers[moltype + "_rev"] = ref_rev_e * Ha_to_eV

        for method, calc in methods.items():
            ts_e = calc["ts"]["output"]["final_energy"]
            rct_e = calc["rct"]["output"]["final_energy"]
            pro_e = calc["pro"]["output"]["final_energy"]
            fwd_e = ts_e - rct_e
            rev_e = ts_e - pro_e

            mae[moltype + "_fwd"][method] = abs(fwd_e - ref_fwd_e) * Ha_to_eV
            mae[moltype + "_rev"][method] = abs(rev_e - ref_rev_e) * Ha_to_eV
            #
            # mae[moltype + "_fwd"][method] = abs(fwd_e - ref_fwd_e)
            # mae[moltype + "_rev"][method] = abs(rev_e - ref_rev_e)

            target_barriers[moltype + "_fwd"][method] = abs(fwd_e) * Ha_to_eV
            target_barriers[moltype + "_rev"][method] = abs(rev_e) * Ha_to_eV

    mae["average"] = {}
    for method in mae["carbonate_fwd"].keys():
        all_mae = []
        for moltype in mae.keys():
            if moltype != "average":
                if method in mae[moltype].keys():
                    all_mae.append(mae[moltype][method])

        avg_mae = np.sum(all_mae) / len(all_mae)
        mae["average"][method] = avg_mae
    return mae


def get_elec_mse(reference, target):
    """

    :param reference:
    :param target:
    :return: MAE: Mean Absolute Error for electronic energies of barrier
    heights as a dict of {"molecule":{"functional": mae, . . . }, . . .,
    "average": {"functional": average_mae, . . .}}
    """
    mae = {}
    ref_barriers = {}
    target_barriers = {}
    for moltype, methods in target.items():

        mae[moltype + "_fwd"] = {}
        mae[moltype + "_rev"] = {}

        target_barriers[moltype + "_fwd"] = {}
        target_barriers[moltype + "_rev"] = {}

        ref_ts_e = reference[moltype]["ts"]["total_corrected"]
        ref_rct_e = reference[moltype]["rct"]["total_corrected"]
        ref_pro_e = reference[moltype]["pro"]["total_corrected"]
        ref_fwd_e = ref_ts_e - ref_rct_e
        ref_rev_e = ref_ts_e - ref_pro_e

        ref_barriers[moltype + "_fwd"] = ref_fwd_e * Ha_to_eV
        ref_barriers[moltype + "_rev"] = ref_rev_e * Ha_to_eV

        for method, calc in methods.items():
            ts_e = calc["ts"]["output"]["final_energy"]
            rct_e = calc["rct"]["output"]["final_energy"]
            pro_e = calc["pro"]["output"]["final_energy"]
            fwd_e = ts_e - rct_e
            rev_e = ts_e - pro_e

            mae[moltype + "_fwd"][method] = (fwd_e - ref_fwd_e) * Ha_to_eV
            mae[moltype + "_rev"][method] = (rev_e - ref_rev_e) * Ha_to_eV

            target_barriers[moltype + "_fwd"][method] = (fwd_e) * Ha_to_eV
            target_barriers[moltype + "_rev"][method] = (rev_e) * Ha_to_eV

    mae["average"] = {}
    for method in mae["carbonate_fwd"].keys():
        all_mae = []
        for moltype in mae.keys():
            if moltype != "average":
                if method in mae[moltype].keys():
                    all_mae.append(mae[moltype][method])

        avg_mae = np.sum(all_mae) / len(all_mae)
        mae["average"][method] = avg_mae
    return mae


def get_cost(data):
    """

    :param reference:
    :param target:
    :return: cost: cputime dict of {"molecule":{"functional": mae, . . . },
    . . ., "average": {"functional": average_mae, . . .}}
    """
    cost = {}
    for moltype, methods in data.items():

        cost[moltype + "_ts"] = {}
        cost[moltype + "_rct"] = {}
        cost[moltype + "_pro"] = {}

        for method, calc in methods.items():
            ts = calc["ts"]["cputime"]
            rct = calc["rct"]["cputime"]
            pro = calc["pro"]["cputime"]

            cost[moltype + "_ts"][method] = ts
            cost[moltype + "_rct"][method] = rct
            cost[moltype + "_pro"][method] = pro

    cost["average"] = {}
    for method in cost["carbonate_ts"].keys():
        all_mae = []
        for moltype in cost.keys():
            if moltype != "average":
                if method in cost[moltype].keys():
                    all_mae.append(cost[moltype][method])

        avg_mae = np.sum(all_mae) / len(all_mae)
        cost["average"][method] = avg_mae
    return cost


def get_thermochemistry_maes(reference, target):
    """

    :param reference:
    :param target:
    :return: MAE: Mean Absolute Error for electronic energies of barrier
    heights as a dict of {"molecule":{"functional": mae, . . . }, . . .,
    "average": {"functional": average_mae, . . .}}
    """
    temp = 298
    g_mae = {}
    s_mae = {}
    h_mae = {}
    for moltype, methods in target.items():

        g_mae[moltype + "_fwd"] = {}
        s_mae[moltype + "_fwd"] = {}
        h_mae[moltype + "_fwd"] = {}
        g_mae[moltype + "_rev"] = {}
        s_mae[moltype + "_rev"] = {}
        h_mae[moltype + "_rev"] = {}

        ref_ts_e = reference[moltype]["ts"]["total_corrected"]
        ref_ts_s = reference[moltype]["ts"]["entropy"] / 1000 * kcal_to_Ha
        ref_ts_h = reference[moltype]["ts"]["enthalpy"] * kcal_to_Ha
        ref_ts_g = ref_ts_e + ref_ts_h  - ref_ts_s * temp

        ref_rct_e = reference[moltype]["rct"]["total_corrected"]
        ref_rct_s = reference[moltype]["rct"]["entropy"] / 1000 * kcal_to_Ha
        ref_rct_h = reference[moltype]["rct"]["enthalpy"] * kcal_to_Ha
        ref_rct_g = ref_rct_e + ref_rct_h - ref_rct_s * temp

        ref_pro_e = reference[moltype]["pro"]["total_corrected"]
        ref_pro_s = reference[moltype]["pro"]["entropy"] / 1000 * kcal_to_Ha
        ref_pro_h = reference[moltype]["pro"]["enthalpy"] * kcal_to_Ha
        ref_pro_g = ref_pro_e + ref_pro_h - ref_pro_s * temp

        ref_fwd_g = ref_ts_g - ref_rct_g
        ref_fwd_s = ref_ts_s - ref_rct_s
        ref_fwd_h = ref_ts_h - ref_rct_h

        ref_rev_g = ref_ts_g - ref_pro_g
        ref_rev_s = ref_ts_s - ref_pro_s
        ref_rev_h = ref_ts_h - ref_pro_h

        for method, calc in methods.items():
            ts_e = calc["ts"]["output"]["final_energy"]
            ts_h = calc["ts"]["output"]["enthalpy"] * kcal_to_Ha
            ts_s = calc["ts"]["output"]["entropy"] / 1000 * kcal_to_Ha
            ts_g = ts_e + ts_h - temp * ts_s

            rct_e = calc["rct"]["output"]["final_energy"]
            rct_h = calc["rct"]["output"]["enthalpy"] * kcal_to_Ha
            rct_s = calc["rct"]["output"]["entropy"] / 1000 * kcal_to_Ha
            rct_g = rct_e + rct_h - temp * rct_s

            pro_e = calc["pro"]["output"]["final_energy"]
            pro_h = calc["pro"]["output"]["enthalpy"] * kcal_to_Ha
            pro_s = calc["pro"]["output"]["entropy"] / 1000 * kcal_to_Ha
            pro_g = pro_e + pro_h - temp * pro_s

            fwd_g = ts_g - rct_g
            fwd_h = ts_h - rct_h
            fwd_s = ts_s - rct_s

            rev_g = ts_g - pro_g
            rev_h = ts_h - pro_h
            rev_s = ts_s - pro_s

            g_mae[moltype + "_fwd"][method] = abs(ref_fwd_g - fwd_g) * Ha_to_eV
            s_mae[moltype + "_fwd"][method] = abs(ref_fwd_s - fwd_s) * Ha_to_eV
            h_mae[moltype + "_fwd"][method] = abs(ref_fwd_h - fwd_h) * Ha_to_eV

            g_mae[moltype + "_rev"][method] = abs(ref_rev_g - rev_g) * Ha_to_eV
            s_mae[moltype + "_rev"][method] = abs(ref_rev_s - rev_s) * Ha_to_eV
            h_mae[moltype + "_rev"][method] = abs(ref_rev_h - rev_h) * Ha_to_eV

    g_mae["average"] = {}
    h_mae["average"] = {}
    s_mae["average"] = {}
    for method in g_mae["carbonate_fwd"].keys():
        all_g_mae = []
        all_h_mae = []
        all_s_mae = []
        for moltype in g_mae.keys():
            if moltype != "average":
                if method in g_mae[moltype].keys():
                    all_g_mae.append(g_mae[moltype][method])
                    all_h_mae.append(h_mae[moltype][method])
                    all_s_mae.append(s_mae[moltype][method])

        avg_g_mae = np.sum(all_g_mae) / len(all_g_mae)
        avg_h_mae = np.sum(all_h_mae) / len(all_h_mae)
        avg_s_mae = np.sum(all_s_mae) / len(all_s_mae)
        g_mae["average"][method] = avg_g_mae
        h_mae["average"][method] = avg_h_mae
        s_mae["average"][method] = avg_s_mae
    return (g_mae, h_mae, s_mae)


def get_thermochemistry_rmaes(reference, target):
    """

    :param reference:
    :param target:
    :return: MAE: Mean Absolute Error for electronic energies of barrier
    heights as a dict of {"molecule":{"functional": mae, . . . }, . . .,
    "average": {"functional": average_mae, . . .}}
    """
    temp = 298
    g_mae = {}
    s_mae = {}
    h_mae = {}
    for moltype, methods in target.items():

        g_mae[moltype + "_fwd"] = {}
        s_mae[moltype + "_fwd"] = {}
        h_mae[moltype + "_fwd"] = {}
        g_mae[moltype + "_rev"] = {}
        s_mae[moltype + "_rev"] = {}
        h_mae[moltype + "_rev"] = {}

        ref_ts_e = reference[moltype]["ts"]["total_corrected"]
        ref_ts_s = reference[moltype]["ts"]["entropy"] / 1000 * kcal_to_Ha
        ref_ts_h = reference[moltype]["ts"]["enthalpy"] * kcal_to_Ha
        ref_ts_g = ref_ts_h + ref_ts_e - ref_ts_s * temp

        ref_rct_e = reference[moltype]["rct"]["total_corrected"]
        ref_rct_s = reference[moltype]["rct"]["entropy"] / 1000 * kcal_to_Ha
        ref_rct_h = reference[moltype]["rct"]["enthalpy"] * kcal_to_Ha
        ref_rct_g = ref_rct_e - ref_rct_s * temp + ref_rct_h

        ref_pro_e = reference[moltype]["pro"]["total_corrected"]
        ref_pro_s = reference[moltype]["pro"]["entropy"] / 1000 * kcal_to_Ha
        ref_pro_h = reference[moltype]["pro"]["enthalpy"] * kcal_to_Ha
        ref_pro_g = ref_pro_e - ref_pro_s * temp + ref_pro_h

        ref_fwd_g = ref_ts_g - ref_rct_g
        ref_fwd_s = ref_ts_s - ref_rct_s
        ref_fwd_h = ref_ts_h - ref_rct_h

        ref_rev_g = ref_ts_g - ref_pro_g
        ref_rev_s = ref_ts_s - ref_pro_s
        ref_rev_h = ref_ts_h - ref_pro_h

        for method, calc in methods.items():
            ts_e = calc["ts"]["output"]["final_energy"]
            ts_h = calc["ts"]["output"]["enthalpy"] * kcal_to_Ha
            ts_s = calc["ts"]["output"]["entropy"] / 1000 * kcal_to_Ha
            ts_g = ts_e + ts_h - temp * ts_s

            rct_e = calc["rct"]["output"]["final_energy"]
            rct_h = calc["rct"]["output"]["enthalpy"] * kcal_to_Ha
            rct_s = calc["rct"]["output"]["entropy"] / 1000 * kcal_to_Ha
            rct_g = rct_e + rct_h - temp * rct_s

            pro_e = calc["pro"]["output"]["final_energy"]
            pro_h = calc["pro"]["output"]["enthalpy"] * kcal_to_Ha
            pro_s = calc["pro"]["output"]["entropy"] / 1000 * kcal_to_Ha
            pro_g = pro_e + pro_h - temp * pro_s

            fwd_g = ts_g - rct_g
            fwd_h = ts_h - rct_h
            fwd_s = ts_s - rct_s

            rev_g = ts_g - pro_g
            rev_h = ts_h - pro_h
            rev_s = ts_s - pro_s

            g_mae[moltype + "_fwd"][method] = abs(ref_fwd_g - fwd_g) / ref_fwd_g
            s_mae[moltype + "_fwd"][method] = abs(ref_fwd_s - fwd_s) / ref_fwd_s
            h_mae[moltype + "_fwd"][method] = abs(ref_fwd_h - fwd_h) / ref_fwd_h

            g_mae[moltype + "_rev"][method] = abs(ref_rev_g - rev_g) / ref_rev_g
            s_mae[moltype + "_rev"][method] = abs(ref_rev_s - rev_s) / ref_rev_s
            h_mae[moltype + "_rev"][method] = abs(ref_rev_h - rev_h) / ref_rev_h

    g_mae["average"] = {}
    for method in g_mae["carbonate_fwd"].keys():
        all_mae = []
        for moltype in g_mae.keys():
            if moltype != "average":
                if method in g_mae[moltype].keys():
                    all_mae.append(g_mae[moltype][method])

        avg_mae = np.sum(all_mae) / len(all_mae)
        g_mae["average"][method] = avg_mae
    return g_mae


def get_thermochemistry_mses(reference, target):
    """

    :param reference:
    :param target:
    :return: MAE: Mean Absolute Error for electronic energies of barrier
    heights as a dict of {"molecule":{"functional": mae, . . . }, . . .,
    "average": {"functional": average_mae, . . .}}
    """
    temp = 298
    g_mae = {}
    s_mae = {}
    h_mae = {}
    for moltype, methods in target.items():

        g_mae[moltype + "_fwd"] = {}
        s_mae[moltype + "_fwd"] = {}
        h_mae[moltype + "_fwd"] = {}
        g_mae[moltype + "_rev"] = {}
        s_mae[moltype + "_rev"] = {}
        h_mae[moltype + "_rev"] = {}

        ref_ts_e = reference[moltype]["ts"]["total_corrected"]
        ref_ts_s = reference[moltype]["ts"]["entropy"] / 1000 * kcal_to_Ha
        ref_ts_h = reference[moltype]["ts"]["enthalpy"] * kcal_to_Ha
        ref_ts_g = ref_ts_h + ref_ts_e - ref_ts_s * temp

        ref_rct_e = reference[moltype]["rct"]["total_corrected"]
        ref_rct_s = reference[moltype]["rct"]["entropy"] / 1000 * kcal_to_Ha
        ref_rct_h = reference[moltype]["rct"]["enthalpy"] * kcal_to_Ha
        ref_rct_g = ref_rct_e - ref_rct_s * temp + ref_rct_h

        ref_pro_e = reference[moltype]["pro"]["total_corrected"]
        ref_pro_s = reference[moltype]["pro"]["entropy"] / 1000 * kcal_to_Ha
        ref_pro_h = reference[moltype]["pro"]["enthalpy"] * kcal_to_Ha
        ref_pro_g = ref_pro_e - ref_pro_s * temp + ref_pro_h

        ref_fwd_g = ref_ts_g - ref_rct_g
        ref_fwd_s = ref_ts_s - ref_rct_s
        ref_fwd_h = ref_ts_h - ref_rct_h

        ref_rev_g = ref_ts_g - ref_pro_g
        ref_rev_s = ref_ts_s - ref_pro_s
        ref_rev_h = ref_ts_h - ref_pro_h

        for method, calc in methods.items():
            ts_e = calc["ts"]["output"]["final_energy"]
            ts_h = calc["ts"]["output"]["enthalpy"] * kcal_to_Ha
            ts_s = calc["ts"]["output"]["entropy"] / 1000 * kcal_to_Ha
            ts_g = ts_e + ts_h - temp * ts_s

            rct_e = calc["rct"]["output"]["final_energy"]
            rct_h = calc["rct"]["output"]["enthalpy"] * kcal_to_Ha
            rct_s = calc["rct"]["output"]["entropy"] / 1000 * kcal_to_Ha
            rct_g = rct_e + rct_h - temp * rct_s

            pro_e = calc["pro"]["output"]["final_energy"]
            pro_h = calc["pro"]["output"]["enthalpy"] * kcal_to_Ha
            pro_s = calc["pro"]["output"]["entropy"] / 1000 * kcal_to_Ha
            pro_g = pro_e + pro_h - temp * pro_s

            fwd_g = ts_g - rct_g
            fwd_h = ts_h - rct_h
            fwd_s = ts_s - rct_s

            rev_g = ts_g - pro_g
            rev_h = ts_h - pro_h
            rev_s = ts_s - pro_s

            g_mae[moltype + "_fwd"][method] = (ref_fwd_g - fwd_g) * Ha_to_eV
            s_mae[moltype + "_fwd"][method] = (ref_fwd_s - fwd_s) * Ha_to_eV
            h_mae[moltype + "_fwd"][method] = (ref_fwd_h - fwd_h) * Ha_to_eV

            g_mae[moltype + "_rev"][method] = (ref_rev_g - rev_g) * Ha_to_eV
            s_mae[moltype + "_rev"][method] = (ref_rev_s - rev_s) * Ha_to_eV
            h_mae[moltype + "_rev"][method] = (ref_rev_h - rev_h) * Ha_to_eV

    g_mae["average"] = {}
    h_mae["average"] = {}
    s_mae["average"] = {}
    for method in g_mae["carbonate_fwd"].keys():
        all_g_mae = []
        all_h_mae = []
        all_s_mae = []
        for moltype in g_mae.keys():
            if moltype != "average":
                if method in g_mae[moltype].keys():
                    all_g_mae.append(g_mae[moltype][method])
                    all_h_mae.append(h_mae[moltype][method])
                    all_s_mae.append(s_mae[moltype][method])

        avg_g_mae = np.sum(all_g_mae) / len(all_g_mae)
        avg_h_mae = np.sum(all_h_mae) / len(all_h_mae)
        avg_s_mae = np.sum(all_s_mae) / len(all_s_mae)
        g_mae["average"][method] = avg_g_mae
        h_mae["average"][method] = avg_h_mae
        s_mae["average"][method] = avg_s_mae
    return (g_mae, h_mae, s_mae)


def get_rmsds(reference, target):
    rmsds = {}
    for moltype, methods in target.items():

        rmsds[moltype + "_ts"] = {}
        rmsds[moltype + "_rct"] = {}
        rmsds[moltype + "_pro"] = {}

        ref_ts_mol = reference[moltype]["ts"]["molecule"]
        ref_rct_mol = reference[moltype]["rct"]["molecule"]
        ref_pro_mol = reference[moltype]["pro"]["molecule"]
        tsmm = KabschMatcher(ref_ts_mol)
        rctmm = KabschMatcher(ref_rct_mol)
        promm = KabschMatcher(ref_pro_mol)
        # mm = MoleculeMatcher()
        for method, calc in methods.items():
            ts_mol = Molecule.from_dict(calc["ts"]["output"]["optimized_molecule"])
            rct_mol = Molecule.from_dict(calc["rct"]["output"]["optimized_molecule"])
            pro_mol = Molecule.from_dict(calc["pro"]["output"]["optimized_molecule"])

            # ts_rmsd = mm.get_rmsd(ts_mol, ref_ts_mol)
            # rct_rmsd = mm.get_rmsd(rct_mol, ref_rct_mol)
            # pro_rmsd = mm.get_rmsd(pro_mol, ref_pro_mol)
            p, ts_rmsd = tsmm.fit(ts_mol)
            p, rct_rmsd = rctmm.fit(rct_mol)
            p, pro_rmsd = promm.fit(pro_mol)
            rmsds[moltype + "_ts"][method] = ts_rmsd
            rmsds[moltype + "_rct"][method] = rct_rmsd
            rmsds[moltype + "_pro"][method] = pro_rmsd
    rmsds["average"] = {}
    for method in rmsds["carbonate_ts"].keys():
        all_mae = []
        for moltype in rmsds.keys():
            if moltype != "average":
                if method in rmsds[moltype].keys():
                    all_mae.append(rmsds[moltype][method])

        avg_mae = np.sum(all_mae) / len(all_mae)
        rmsds["average"][method] = avg_mae
    return rmsds


def get_svpd_with_sp_corrections(data):
    spcorr_data = copy.deepcopy(data)
    for moltype, methods in data.items():
        for method, calc in methods.items():
            if "ts_sp" in calc.keys():
                ts_e = calc["ts_sp"]["output"]["final_energy"]
                spcorr_data[moltype][method]["ts"]["output"]["final_energy"] = ts_e
                rct_e = calc["rct_sp"]["output"]["final_energy"]
                spcorr_data[moltype][method]["rct"]["output"]["final_energy"] = rct_e
                pro_e = calc["pro_sp"]["output"]["final_energy"]
                spcorr_data[moltype][method]["pro"]["output"]["final_energy"] = pro_e
            else:
                print("Missing single point for {},{}".format(moltype, method))
                spcorr_data[moltype].pop(method)

    return spcorr_data


mp2_ref = loadfn(os.path.join(data_dir, "mp2opt_svpd_reference_data.json"))
# fix amide name
mp2_ref["amide"] = mp2_ref.pop("amide22")
with open(os.path.join(data_dir, "opt_data_dump_tzvppd.json")) as json_file:
    tzvpd_data = json.load(json_file)
with open(os.path.join(data_dir, "opt_data_dump_svpd.json")) as json_file:
    svpd_data = json.load(json_file)

elec_mae = get_elec_mae(reference=mp2_ref, target=tzvpd_data)
plot_heatmap(elec_mae, title='E_MAE_TZVPPD')
print("Elec MAE TZVPD")
dumpfn(elec_mae["average"], fn="../data/elec_mae_tzvpd_opts.json")
g_mae, h_mae, s_mae = get_thermochemistry_maes(reference=mp2_ref, target=tzvpd_data)
# print('G, H, S MAE TZVPD')
dumpfn(g_mae["average"], fn="../data/g_mae_tzvpd_opts.json")
dumpfn(h_mae["average"], fn="../data/h_mae_tzvpd_opts.json")
dumpfn(s_mae["average"], fn="../data/s_mae_tzvpd_opts.json")

plot_heatmap(g_mae, title='G_MAE_TZVPPD')

g_rmae = get_thermochemistry_rmaes(reference=mp2_ref, target=tzvpd_data)
# plot_heatmap(g_rmae, title='G_RMAE_TZVPPD')

g_mae_svpd, h_mae_svpd, s_mae_svpd = get_thermochemistry_maes(
    reference=mp2_ref, target=svpd_data
)
# plot_heatmap(g_mae_svpd, title='G_MAE_SVPD')
dumpfn(g_mae_svpd["average"], fn="../data/g_mae_svpd_opts.json")
dumpfn(h_mae_svpd["average"], fn="../data/h_mae_svpd_opts.json")
dumpfn(s_mae_svpd["average"], fn="../data/s_mae_svpd_opts.json")
# avg_g_mae_svpd

g_rmae_svpd = get_thermochemistry_rmaes(reference=mp2_ref, target=svpd_data)
# plot_heatmap(g_rmae_svpd, title='G_RMAE_SVPD')

tzvppd_cost = get_cost(data=tzvpd_data)
# plot_heatmap(tzvppd_cost)

rmsds = get_rmsds(reference=mp2_ref, target=tzvpd_data)
# plot_heatmap(rmsds,zmax=None)

spcorr_svpd = get_svpd_with_sp_corrections(data=svpd_data)

elec_mse = get_elec_mse(reference=mp2_ref, target=spcorr_svpd)
# plot_heatmap(elec_mse, title='Elec_E_MSE')


g_mae_corrsvpd, h_mae_corrsvpd, s_mae_corrsvpd = get_thermochemistry_maes(
    reference=mp2_ref, target=spcorr_svpd
)
dumpfn(g_mae_corrsvpd["average"], fn="../data/g_mae_corrsvpd_opts.json")
dumpfn(h_mae_corrsvpd["average"], fn="../data/h_mae_corrsvpd_opts.json")
dumpfn(s_mae_corrsvpd["average"], fn="../data/s_mae_corrsvpd_opts.json")

#
g_mse_corrsvpd, h_mse_corrsvpd, s_mse_corrsvpd = get_thermochemistry_mses(
    reference=mp2_ref, target=spcorr_svpd
)

g_rmae_corrsvpd = get_thermochemistry_rmaes(reference=mp2_ref, target=spcorr_svpd)
# plot_heatmap(g_rmae_corrsvpd, title='G_RMAE_corrSVPD')

exit()
