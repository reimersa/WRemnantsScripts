import argparse
import os
import pickle
import re

import lz4.frame

from utilities import common
from wremnants import theory_corrections
from wums import boostHistHelpers as hh
from wums import logging, output_tools, plot_tools
from scripts.corrections.make_theory_corr import read_corr
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import mplhep
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
import hist


def corr_name(corrf):
    if not corrf.endswith(".pkl.lz4"):
        raise ValueError(f"File {corrf} is not a lz4 compressed pickle file")

    match = re.match(r"(.*)Corr[Z|W]\.pkl\.lz4", os.path.basename(corrf))
    return match[1]



# Define the exponential model: A * exp(-ax)
def exp_func(x, A, a):
    return A * np.exp(a*x)


def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections"
    
    filenames_resum = [
        # ("SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl", "N$^{4+0}$LL (for $\\alpha_{s}$)"),

        # ("SCETlib/ct18z_nplambda_n4+0ll_arnetest/inclusive_Z_CT18Z_nplambda_N4+0LL_combined.pkl", "N$^{4+0}$LL (made by Arne)"),
        # ("SCETlib/ct18z_nplambda_n4+0ll_arnetest_asn4ll/inclusive_Z_CT18Z_nplambda_N4+0LL_combined.pkl", "N$^{4+0}$LL (as@N$^{4}$LL, made by Arne)"),
        # ("SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_N4+0LL_combined.pkl", "N$^{4+0}$LL (made by Kenneth)"),

        # ("SCETlib/ct18z_nplambda_n3+1ll/inclusive_Z_CT18Z_nplambda_N3+1LL_combined.pkl", "N$^{3+1}$LL"),
        # ("SCETlib/ct18z_nplambda_n3+1ll_arnetest/inclusive_Z_CT18Z_nplambda_N3+1LL_combined.pkl", "N$^{3+1}$LL (made by Arne)"),
        
        # ("SCETlib/ct18z_nplambda_n3+0ll/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl", "N$^{3+0}$LL (made by Kenneth)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_arnetest/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl", "N$^{3+0}$LL (made by Arne)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps/inclusive_Z_CT18Z_N3+0LL_zeros.pkl", "N$^{3+0}$LL (zeros)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps/inclusive_Z_CT18Z_N3+0LL_frank.pkl", "N$^{3+0}$LL (frank)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps/inclusive_Z_CT18Z_N3+0LL_reproducesold.pkl", "N$^{3+0}$LL (reproduces)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps/inclusive_Z_CT18Z_N3+0LL_testing.pkl", "N$^{3+0}$LL (testing)"),

        ("SCETlib/ct18z_nplambda_n3+0ll/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl", "N$^{3+0}$LL (old NP model)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2.pkl", "New NP, repr. old"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_lambdainfnu_2p0.pkl", "+ $\\lambda_{\\infty}^{\\nu}$ = +2.0"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_lambda4nu_0p01.pkl", "repr. + $\\lambda_{4}^{\\nu}$ = 0.01"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_b0overbmaxnu_1p0.pkl", "repr. + b0/bmax = 1.0"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_lambdainfnu_m0p2.pkl", "(orig. + $\\lambda_{\\infty}^{\\nu}$ = -0.2)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2.pkl", "(orig. + NP model tanh_2)"),

        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_m0p2.pkl", "new NP reproducing old NP"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2_0p25_combined.pkl", "$\\lambda_2=0.25$ (standalone check)"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2nu_0p1_combined.pkl", "+ $\\lambda_2^{\\nu}$=0.1"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2nu_0p1_lambda4nu_0p01_combined.pkl", "+ $\\lambda_4^{\\nu}$=0.01"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2nu_0p1_lambda4nu_0p01_b0overbmaxnu_1p0_combined.pkl", "+ b$_{0}$/b$_{max}$=1.0"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2nu_0p1_lambda4nu_0p01_b0overbmaxnu_1p0_lambda2_0p25_combined.pkl", "+ $\\lambda_2=0.25$"),
        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps_stepbystep/inclusive_Z_CT18Z_N3+0LL_novars_npmodel_tanh2_lambdainfnu_2p0_lambda2nu_0p1_lambda4nu_0p01_b0overbmaxnu_1p0_lambda2_0p25_deltalambda2_0p125.pkl", "+ $\\Delta\\lambda_2=0.125$"),


        # ("SCETlib/ct18z_nplambda_n3+0ll_newnps/inclusive_Z_CT18Z_N3+0LL_frank.pkl", "N$^{3+0}$LL (Frank)"),
        ("SCETlib/ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl", "N$^{3+0}$LL (lattice)"),


        ("SCETlib/ct18z_newnps_n3+0ll_reproducesold/inclusive_Z_CT18Z_N3+0LL_olduncertainties_combined.pkl", "N$^{3+0}$LL (new NP reproducing old)"),
        
    ]
    
    filenames_fosing = [
        ("SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl", "FO sing. NNLO"),
        ("SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_n3lo_sing.pkl", "FO sing. N$^{3}$LO"),
    ]

    filenames_fo = [
        ("DYTURBO/nnlo-scetlibmatch/scalevariations/z0/results_z-2d-nnlo-vj-CT18ZNNLO-{scale}-scetlibmatch.txt", "FO NNLO"),
        ("NNLOjet/Z/ZjNNLO/final/ptz", "FO N$^{3}$LO"),
    ]

    hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, fnr[0]))[{"vars" :  0}] for fnr in filenames_resum]
    hists_fosing = [input_tools.read_scetlib_hist(os.path.join(genpath, fns[0]))[{"vars" :  0}] for fns in filenames_fosing]
    hists_fo = []
    for idx, (fnf, _) in enumerate(filenames_fo):
        if fnf.startswith("DYTURBO"):
            hists_fo.append(input_tools.read_dyturbo_vars_hist(os.path.join(genpath, fnf), var_axis=hist.axis.StrCategory(['pdf0', 'kappaFO0.5-kappaf2.', 'kappaFO2.-kappaf0.5', 'kappaf0.5', 'kappaf2.', 'kappaFO0.5', 'kappaFO2.'], name="vars"), axes=['Y', 'qT'], charge=0)[{"vars" :  0}])
        elif fnf.startswith("NNLOjet"):
            y_axes = np.array([-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ])
            hists_fo.append(input_tools.read_nnlojet_pty_hist(os.path.join(genpath, fnf), ybins=y_axes, charge=0)[{"vars" :  0}])
        else:
            raise ValueError(f"invalid FO correction, does not start with DYTURBO or SCETLIB: {fnf}")

    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    os.makedirs(outfoldername, exist_ok=True)


    variables = ["qT", "Y"]
    # variables = ["qT"]

    for v in variables:
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[hist_resum.project(v) for hist_resum in hists_resum],
            labels=[fnr[1] for fnr in filenames_resum],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(filenames_resum)],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="Resumm. prediction",
            rlabel=f"x/{filenames_resum[0][1]}",
            rrange=[0.80, 1.10] if v=="qT" else [0.9995, 1.0002],
            nlegcols=1,
            # xlim=None, 
            xlim=[0., 100.] if v == "qT" else [-5., 5.], 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_comparison_resum_{v}.pdf")

        if v=="Y":
            continue

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[hist_fosing.project(v) for hist_fosing in hists_fosing],
            labels=[fns[1] for fns in filenames_fosing],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(filenames_fosing)],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="FO sing. prediction",
            rlabel=f"x/{filenames_fosing[0][1]}",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_comparison_fosingular_{v}.pdf")

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[hist_fo.project(v) for hist_fo in hists_fo],
            labels=[fnf[1] for fnf in filenames_fo],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(filenames_fo)],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="FO. prediction",
            rlabel=f"x/{filenames_fo[0][1]}",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_comparison_fo_{v}.pdf")

        if len(filenames_fosing) == len(filenames_fo):
            hists_nonsing = []
            for idxh, h in enumerate([hfo.project(v).copy() for hfo in hists_fo]):
                h.view(flow=False)[...] = hists_fo[idxh].project("qT").view(flow=False) - hists_fosing[idxh].project("qT").view(flow=False)
                hists_nonsing.append(h)

            # fo_nonsing = var_fo.project("qT").copy()
            #fo_nonsing.view(flow=False)[...] = var_fo.project("qT").view(flow=False) - var_fosing.project("qT").view(flow=False)


            fig = plot_tools.makePlotWithRatioToRef(
            hists=[hist_nonsing.project(v) for hist_nonsing in hists_nonsing],
            labels=[fns[1].replace("sing", "nonsing") for fns in filenames_fosing],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(filenames_fosing)],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="FO nonsing. prediction",
            rlabel=f"x/{filenames_fosing[0][1]}",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_comparison_fononsingular_{v}.pdf")

if __name__ == "__main__":
    main()