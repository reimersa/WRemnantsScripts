import os
import re


from wums import plot_tools
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors


def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections/SCETlib"
    
    # filename_resum = "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl"

    variations_and_values_per_par = {
        "\\Lambda_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            "pdf0": 0.25,
            "lambda20.0": 0.0,
            "lambda20.5": 0.5,
            "lambda2-0.25": -0.25,
            "lambda20.75": 0.75,
            "lambda2-0.5": -0.5,
            "lambda21.0": 1.0,
        }),
        "\\Lambda_{4}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.0625,
            # "lambda40.0": 0.0,
            # "lambda40.001": 0.001,
            # "lambda40.125": 0.125,

            # for lattice
            "pdf0": 0.06,
            'lambda40.0': 0.0, 
            'lambda40.12': 0.12, 
            'lambda40.04': 0.04, 
            'lambda40.08': 0.08
        }),
        "\\Delta\\Lambda_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.125,
            # "delta_lambda20.105": 0.105,
            # "delta_lambda20.145": 0.145,
            
            # for lattice
            "pdf0": 0.125,
            'delta_lambda20.105': 0.105, 
            'delta_lambda20.145': 0.145, 
            'delta_lambda20.085': 0.085, 
            'delta_lambda20.165': 0.165, 
            'delta_lambda20.065': 0.065, 
            'delta_lambda20.185': 0.185
        }),
        "\\lambda^{\\nu}_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.1,
            # "lambda2_nu0.0": 0.0,
            # "lambda2_nu0.2": 0.2,
            # "lambda2_nu-0.1": -0.1,
            # "lambda2_nu0.3": 0.3,
            # "lambda2_nu-0.2": -0.2,
            # "lambda2_nu0.4": 0.4,
            # "lambda2_nu-0.3": -0.3,
            # "lambda2_nu0.5": 0.5,
            # "lambda2_nu-0.4": -0.4,
            # "lambda2_nu0.6": 0.6,

            # for lattice
            "pdf0": 0.087,
            'lambda2_nu0.0538': 0.0538, 
            'lambda2_nu0.1202': 0.1202

        }),
        "\\lambda^{\\nu}_{4}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl",
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.01,
            # "lambda4_nu0.0": 0.0,
            # "lambda4_nu0.0001": 0.0001,
            # "lambda4_nu0.02": 0.02,
            # "lambda4_nu0.03": 0.03,
            # "lambda4_nu0.006": 0.006,

            # for lattice
            "pdf0": 0.0074,
            'lambda4_nu0.0008': 0.0008, 
            'lambda4_nu0.014': 0.014
        }),
    }

    qts = [0, 5, 10, 15, 20, 40]
    # qts = [0]


    for par in variations_and_values_per_par:
        print(input_tools.read_scetlib_hist(os.path.join(genpath, variations_and_values_per_par[par][0])))

        hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, variations_and_values_per_par[par][0]))[{"vars" : var}] for var in variations_and_values_per_par[par][1]]
        outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
        os.makedirs(outfoldername, exist_ok=True)

        values_per_qt = {qt: [(list(variations_and_values_per_par[par][1].items())[idx][1], hists_resum[idx].project("qT").values()[hists_resum[idx].project("qT").axes[0].index(qt)]) for idx in range(len(variations_and_values_per_par[par][1]))] for qt in qts}

        plt.figure(figsize=(7,5))

        for key, points in values_per_qt.items():
            # sort by x-value
            points_sorted = sorted(points, key=lambda t: t[0])
            xs, ys = zip(*points_sorted)
            y_ref = next(y for x,y in points_sorted if x == variations_and_values_per_par[par][1]["pdf0"])
            ys_shifted = [y - y_ref for y in ys]
            plt.plot(xs, ys_shifted, marker="o", label=f"$q_{{T}}={key}$ GeV")

        plt.xlabel(f"${par}$", fontsize=20)
        plt.ylabel("Resummed prediction - nominal", fontsize=20)
        plt.tick_params(labelsize=18)
        plt.legend(fontsize=18)
        plt.tight_layout()
        plt.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.15)
        plt.savefig(f"{outfoldername}/correction_variations_vs_parameter_{par.replace("\\", "").replace("_", "").replace("{", "").replace("}", "").replace("^", "")}_resum.pdf")

        plt.xscale("log")
        plt.savefig(f"{outfoldername}/correction_variations_vs_parameter_log_{par.replace("\\", "").replace("_", "").replace("{", "").replace("}", "").replace("^", "")}_resum.pdf")

        

if __name__ == "__main__":
    main()