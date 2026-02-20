import os
import re


from wums import plot_tools
import matplotlib.pyplot as plt # type: ignore
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors # type: ignore
import itertools

from wums import boostHistHelpers as hh

def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections"
    filenames_to_compare_from = [
        ("SCETlib/ct18z_newnps_n3+0ll_lattice_pdfvars/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_withnpvars_pdfas_pdf_CT18ZNNLO_as_0118_combined.pkl", "as = 0.118"),
        ("SCETlib/ct18z_newnps_n3+0ll_lattice_pdfvars/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_withnpvars_pdfas_pdf_CT18ZNNLO_as_0116_combined.pkl", "as = 0.116"),
        ("SCETlib/ct18z_newnps_n3+0ll_lattice_pdfvars/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_withnpvars_pdfas_pdf_CT18ZNNLO_as_0120_combined.pkl", "as = 0.120"),
    ]
    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone/variations_per_alphas"
    

    # will make one plot per item in the outer list, comparing the files in the list in each plot
    variations_to_plot = {

        # pdfas vars WITH varied NP parameters (= around different central values)
        "Eigvar1": [
            ("pdf0", "nominal"),
            ("lambda2_nu0.0696-lambda4_nu0.0122-lambda_inf_nu1.1Ext", "eigvar 1 up"),
            ("lambda2_nu0.1044-lambda4_nu0.0026-lambda_inf_nu2.1Ext", "eigvar 1 down")
        ],
        "Eigvar2": [
            ("pdf0", "nominal"),
            ("lambda2_nu0.1153-lambda4_nu0.0032-lambda_inf_nu1.6Ext", "eigvar 2 up"),
            ("lambda2_nu0.0587-lambda4_nu0.0116-lambda_inf_nu1.6Ext", "eigvar 2 down"),
        ],
        "Eigvar3": [
            ("pdf0", "nominal"),
            ("lambda2_nu0.0873-lambda4_nu0.0092", "eigvar 3 up"),
            ("lambda2_nu0.0867-lambda4_nu0.0056", "eigvar 3 down"),
        ],
        "Lambda2": [
            ("pdf0", "nominal"),
            ("lambda20.0", "$\\Lambda_{2}^{new} = 0.00$"),
            ("lambda20.5", "$\\Lambda_{2}^{new} = 0.50$"),
        ],
        "Lambda4": [
            ("pdf0", "nominal"),
            ("lambda40.01", "$\\Lambda_{4}^{new} = 0.01$"),
            ("lambda40.16", "$\\Lambda_{4}^{new} = 0.16$"),
        ],
        "DeltaLambda2": [
            ("pdf0", "nominal"),
            ("delta_lambda20.105", "$\\Delta\\Lambda_{2}^{new} = 0.105$"),
            ("delta_lambda20.145", "$\\Delta\\Lambda_{2}^{new} = 0.145$"),
        ],
        
    }




    # print(input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum)))

    os.makedirs(outfoldername, exist_ok=True)
    for variationname in variations_to_plot.keys():
        hists = [input_tools.read_scetlib_hist(os.path.join(genpath, fn_legas[0]))[{"vars" : var_legvar[0]}] for (fn_legas, var_legvar) in list(itertools.product(filenames_to_compare_from, variations_to_plot[variationname]))]
        denhists = [input_tools.read_scetlib_hist(os.path.join(genpath, filenames_to_compare_from[0][0]))[{"vars" : var_legvar[0]}] for (fn_legas, var_legvar) in list(itertools.product(filenames_to_compare_from, variations_to_plot[variationname]))] # always make ratio to first in list of files to compare (as = 0.118)


        variables = ["qT", "Y"]

        for v in variables:
            
            ratio_hists = [
                hh.divideHists(
                    h.project(v),
                    hden.project(v),
                    cutoff=1e-6,
                    flow=False,
                    rel_unc=True,
                    by_ax_name=False,
                )
                for h, hden in zip(hists, denhists)
            ]



            fig = plot_tools.makePlotWithRatioToRef(
                hists=[hist.project(v) for hist in hists],
                hists_ratio=ratio_hists,
                override_hists_ratio=True,
                labels=[f"{fn_legas[1]}, {var_legvar[1]}" for (fn_legas, var_legvar) in list(itertools.product(filenames_to_compare_from, variations_to_plot[variationname]))],
                colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab20").colors)][:len(list(itertools.product(filenames_to_compare_from, variations_to_plot[variationname])))],
                xlabel=f"{v} (GeV)" if v=="qT" else v, 
                ylabel="Resumm. prediction",
                rlabel=f"x/0.118",
                rrange=[0.95, 1.05] if v=="qT" else [0.9990, 1.0010],
                nlegcols=1,
                xlim=[0., 100.] if v == "qT" else [-5., 5.], 
                binwnorm=1.0, 
                baseline=True, 
                yerr=True,
                yerr_ratio=True,
                ratio_legend=False,
                linewidth=1,
                logy=False,
            )
            fig.savefig(f"{outfoldername}/{variationname}_{v}.pdf")
            print(f"made plot '{outfoldername}/{variationname}_{v}.pdf'")
            plt.close(fig)


        

if __name__ == "__main__":
    main()