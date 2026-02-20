import os
import re


from wums import plot_tools
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors


def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections/SCETlib"
    
    variations_and_values_per_par = {
        "\\Lambda_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            "pdf0": 0.25,
            # "lambda20.0": 0.0,
            # "lambda20.5": 0.5,
            # "lambda2-0.25": -0.25,
            # "lambda20.75": 0.75,
            # "lambda2-0.5": -0.5,
            # "lambda21.0": 1.0,

            "lambda2-0.5": -0.5, 
            "lambda2-0.25": -0.25, 
            "lambda20.0": 0.0, 
            "lambda20.06": 0.06, 
            "lambda20.125": 0.125, 
            "lambda20.375": 0.375, 
            "lambda20.44": 0.44, 
            "lambda20.5": 0.5, 
            "lambda20.75": 0.75, 
            "lambda21.0": 1.0, 
        }),
        "\\Lambda_{4}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.0625,
            # "lambda40.0": 0.0,
            # "lambda40.001": 0.001,
            # "lambda40.125": 0.125,

            # for lattice
            "pdf0": 0.06,
            # 'lambda40.0': 0.0, 
            # 'lambda40.12': 0.12, 
            # 'lambda40.04': 0.04, 
            # 'lambda40.08': 0.08

            "lambda4-0.24": -0.24, 
            "lambda4-0.18": -0.18, 
            "lambda4-0.12": -0.12, 
            "lambda4-0.06": -0.06, 
            "lambda40.00": 0.00, 
            "lambda40.01": 0.01, 
            "lambda40.02": 0.02, 
            "lambda40.04": 0.04, 
            "lambda40.08": 0.08, 
            "lambda40.10": 0.10, 
            "lambda40.12": 0.12, 
            "lambda40.16": 0.16, 
            "lambda40.18": 0.18, 
            "lambda40.24": 0.24, 
            "lambda40.30": 0.30, 
            "lambda40.36": 0.36, 
        }),
        "\\Delta\\Lambda_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # "pdf0": 0.125,
            # "delta_lambda20.105": 0.105,
            # "delta_lambda20.145": 0.145,
            
            # for lattice
            "pdf0": 0.125,
            # 'delta_lambda20.105': 0.105, 
            # 'delta_lambda20.145': 0.145, 
            # 'delta_lambda20.085': 0.085, 
            # 'delta_lambda20.165': 0.165, 
            # 'delta_lambda20.065': 0.065, 
            # 'delta_lambda20.185': 0.185


            "delta_lambda2-0.075": -0.075, 
            "delta_lambda2-0.055": -0.055, 
            "delta_lambda2-0.035": -0.035, 
            "delta_lambda2-0.015": -0.015, 
            "delta_lambda20.000": 0.000, 
            "delta_lambda20.005": 0.005, 
            "delta_lambda20.025": 0.025, 
            "delta_lambda20.045": 0.045, 
            "delta_lambda20.065": 0.065, 
            "delta_lambda20.085": 0.085, 
            "delta_lambda20.105": 0.105, 
            "delta_lambda20.115": 0.115, 
            "delta_lambda20.120": 0.120, 
            "delta_lambda20.130": 0.130, 
            "delta_lambda20.135": 0.135, 
            "delta_lambda20.145": 0.145, 
            "delta_lambda20.165": 0.165, 
            "delta_lambda20.185": 0.185, 
            "delta_lambda20.205": 0.205, 
            "delta_lambda20.225": 0.225, 
            "delta_lambda20.265": 0.265, 
            "delta_lambda20.325": 0.325, 
        }),
        "\\Lambda_{\\infty}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {
            # for lattice
            "pdf0": 1.0,

            "lambda_inf-1.0": -1.0, 
            "lambda_inf0.0": 0.0, 
            "lambda_inf0.5": 0.5, 
            "lambda_inf0.8": 0.8, 
            "lambda_inf0.9": 0.9, 
            "lambda_inf1.1": 1.1, 
            "lambda_inf1.2": 1.2, 
            "lambda_inf1.5": 1.5, 
            "lambda_inf2.0": 2.0, 
            "lambda_inf3.0": 3.0, 
        }),
        "\\lambda^{\\nu}_{2}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", 
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
            # 'lambda2_nu0.0538': 0.0538, 
            # 'lambda2_nu0.1202': 0.1202

            "lambda2_nu0.0100": 0.0100, 
            "lambda2_nu0.0206": 0.0206, 
            "lambda2_nu0.0538": 0.0538, 
            "lambda2_nu0.0000": 0.0000, 
            "lambda2_nu0.1202": 0.1202, 
            "lambda2_nu0.1534": 0.1534, 
            "lambda2_nu0.1866": 0.1866, 
            "lambda2_nu0.2198": 0.2198, 
            "lambda2_nu0.2530": 0.2530, 
            "lambda2_nu0.3360": 0.3360, 
            "lambda2_nu0.4190": 0.4190, 
            "lambda2_nu0.5850": 0.5850, 
            "lambda2_nu0.7510": 0.7510, 
            "lambda2_nu1.0830": 1.0830
        }),
        "\\lambda^{\\nu}_{4}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl",
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
            # 'lambda4_nu0.0008': 0.0008, 
            # 'lambda4_nu0.014': 0.014

            "lambda4_nu0.00575": 0.00575 , 
            "lambda4_nu0.0041":  0.0041 , 
            "lambda4_nu0.0008":  0.0008 , 
            "lambda4_nu0.0000":  0.0000 , 
            "lambda4_nu0.014":   0.014 , 
            "lambda4_nu0.0206":  0.0206 , 
            "lambda4_nu0.0272":  0.0272 , 
            "lambda4_nu0.0338":  0.0338 , 
            "lambda4_nu0.0404":  0.0404 , 
            "lambda4_nu0.0569":  0.0569 , 
            "lambda4_nu0.0734":  0.0734 , 
            "lambda4_nu0.1064":  0.1064 , 
            "lambda4_nu0.1394":  0.1394 , 
            "lambda4_nu0.2054":  0.2054 
        }),
        "\\lambda^{\\nu}_{\\infty}": (
            "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_largervars_higherprecision_combined.pkl", 
            # "ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl",
            # "ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl",
        {

            # for lattice
            "pdf0": 1.6853,

            # "lambda_inf_nu0.1646": 0.1646, 
            "lambda_inf_nu0.6715": 0.6715, 
            "lambda_inf_nu1.1784": 1.1784, 
            "lambda_inf_nu2.1922": 2.1922, 
            "lambda_inf_nu2.6991": 2.6991, 
            "lambda_inf_nu3.2060": 3.2060, 
            "lambda_inf_nu3.7129": 3.7129, 
            "lambda_inf_nu4.2198": 4.2198, 
            "lambda_inf_nu5.48705": 5.48705, 
            "lambda_inf_nu6.7543": 6.7543, 
            "lambda_inf_nu9.2888": 9.2888, 
            "lambda_inf_nu11.8233": 11.8233, 
            "lambda_inf_nu16.8923": 16.8923
        }),
    }

    lattice_upvar_per_par = {
        "\\Lambda_{2}": 0.5, 
        "\\Lambda_{4}": 0.16, 
        "\\Delta\\Lambda_{2}": 0.145, 
        "\\Lambda_{\\infty}": 1.0,
        "\\lambda^{\\nu}_{2}": 0.1202, 
        "\\lambda^{\\nu}_{4}": 0.014, 
        "\\lambda^{\\nu}_{\\infty}": 2.1922, 
    }

    qts = [0, 2, 5, 10, 15, 20]
    # qts = [0]

    modes = ["lin", "logx"]


    for mode in modes:
        for par in variations_and_values_per_par:
            # print(input_tools.read_scetlib_hist(os.path.join(genpath, variations_and_values_per_par[par][0])))

            hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, variations_and_values_per_par[par][0]))[{"vars" : var}] for var in variations_and_values_per_par[par][1]]
            outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
            os.makedirs(outfoldername, exist_ok=True)

            values_per_qt = {qt: [(list(variations_and_values_per_par[par][1].items())[idx][1], hists_resum[idx].project("qT").values()[hists_resum[idx].project("qT").axes[0].index(qt)]) for idx in range(len(variations_and_values_per_par[par][1]))] for qt in qts}

            fig, (ax_main, ax_ratio) = plt.subplots(2, 1, figsize=(7,6), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
            for ax in (ax_main, ax_ratio):
                for spine in ax.spines.values():
                    spine.set_linewidth(1.3)
            for key, points in values_per_qt.items():
                # sort by x-value
                points_sorted = sorted(points, key=lambda t: t[0])
                xs, ys = zip(*points_sorted)
                xs = np.array(xs)
                ys = np.array(ys)
                x_nom = variations_and_values_per_par[par][1]["pdf0"]
                x_up_lattice = lattice_upvar_per_par[par]
                delta_up_lattice = x_up_lattice - x_nom
                y_ref = next(y for x,y in points_sorted if x == x_nom)
                ys_shifted = np.array([(y)/y_ref for y in ys])
                line, = ax_main.plot(xs, ys_shifted, marker="o", markersize=4, linestyle="none", label=f"$q_{{T}}={key}$ GeV")
                color = line.get_color()

                if mode == "lin":
                    xmin_fit = max(0.0, xs.min())
                    xmax_fit = min(x_nom+10*delta_up_lattice, xs.max())
                    mask_fit = (xs >= xmin_fit) & (xs <= xmax_fit)
                    m, b = np.polyfit(xs[mask_fit], ys_shifted[mask_fit], 1)
                    x_fit = np.linspace(xmin_fit, xmax_fit, 100)
                    y_fit = m * x_fit + b
                    y_fit_at_xs = m * xs + b
                    ratio = ys_shifted / y_fit_at_xs
                    ax_main.plot(x_fit, y_fit, linestyle="--", color=color)
                elif mode == "logx":
                    eps = 1e-3
                    xmin_fit = max(eps, xs.min())
                    xmax_fit = min(x_nom + 10 * delta_up_lattice, xs.max())

                    mask_fit = (xs >= xmin_fit) & (xs <= xmax_fit) & (xs > 0) & (ys_shifted > 0)

                    lx = np.log(xs[mask_fit])
                    ly = np.log(ys_shifted[mask_fit])

                    m, b = np.polyfit(lx, ly, 1)  # ln y = m ln x + b
                    x_fit = np.logspace(np.log10(xmin_fit), np.log10(xmax_fit), 200)
                    y_fit = np.exp(b) * x_fit**m

                    # evaluate fit at all xs (for ratio); only valid where xs>0
                    y_fit_at_xs = np.full_like(xs, np.nan, dtype=float)
                    valid = xs > 0
                    y_fit_at_xs[valid] = np.exp(b) * xs[valid]**m

                    ratio = ys_shifted / y_fit_at_xs

                    ax_main.set_xscale("log")
                    ax_ratio.set_xscale("log")
                    ax_main.plot(x_fit, y_fit, linestyle="--", color=color)

                ax_main.axvline(x=x_nom+1*delta_up_lattice, color="black", linestyle="--", linewidth=1.2)
                ax_main.text(x_nom+1*delta_up_lattice, 1.13, "$+1\\sigma$", rotation=90, va="top", ha="right", fontsize=14)
                ax_main.axvline(x=x_nom+10*delta_up_lattice, color="black", linestyle="--", linewidth=1.2)
                ax_main.text(x_nom+10*delta_up_lattice, 1.13, "$+10\\sigma$", rotation=90, va="top", ha="right", fontsize=14)
                ax_main.axvline(x=x_nom, color="black", linestyle="-", linewidth=1.2)
                ax_main.text(x_nom, 1.13, "nominal", rotation=90, va="top", ha="right", fontsize=14)

                ax_ratio.plot(xs, ratio, marker="o", linestyle="none", markersize=4, color=color)
            ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1)

            ax_ratio.set_xlabel(f"${par}$", fontsize=20)
            ax_main.set_ylabel("var / nom", fontsize=20)
            ax_ratio.set_ylabel("points / fit", fontsize=20)
            ax_main.tick_params(labelsize=18)
            ax_ratio.tick_params(labelsize=18)
            ax_main.legend(fontsize=18)
            # fig.tight_layout()
            plt.subplots_adjust(left=0.14, right=0.98, top=0.98, bottom=0.14, hspace=0.12)
            if mode == "logx":
                ax_main.set_xscale("log")
            ax_main.set_ylim(0.80, 1.20)
            ax_ratio.set_ylim(0.95, 1.05)
            plt.savefig(f"{outfoldername}/correction_variations_vs_parameter_{par.replace("\\", "").replace("_", "").replace("{", "").replace("}", "").replace("^", "")}_{mode}_resum.pdf")
            ax_main.set_ylim(0.97, 1.03)
            ax_ratio.set_ylim(0.99, 1.01)
            plt.savefig(f"{outfoldername}/correction_variations_vs_parameter_{par.replace("\\", "").replace("_", "").replace("{", "").replace("}", "").replace("^", "")}_{mode}_resum_zoom.pdf")

            # plt.savefig(f"{outfoldername}/correction_variations_vs_parameter_log_{par.replace("\\", "").replace("_", "").replace("{", "").replace("}", "").replace("^", "")}_resum.pdf")

        

if __name__ == "__main__":
    main()