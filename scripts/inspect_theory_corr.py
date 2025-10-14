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


def corr_name(corrf):
    if not corrf.endswith(".pkl.lz4"):
        raise ValueError(f"File {corrf} is not a lz4 compressed pickle file")

    match = re.match(r"(.*)Corr[Z|W]\.pkl\.lz4", os.path.basename(corrf))
    return match[1]



# Define the exponential model: A * exp(-ax)
def exp_func(x, A, a):
    return A * np.exp(a*x)


def main():


    correctionname1 = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_ProperPdfAsCorrZ.pkl.lz4"
    # correctionname1 = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4"
    # correctionname2 = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfasCorrZ.pkl.lz4"
    correctionname2 = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfasCorrZ.pkl.lz4"
    corrname1=corr_name(correctionname1)
    corrname2=corr_name(correctionname2)
    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    plotoutname_base = f"{corrname1}TESTCORR"



    os.makedirs(outfoldername, exist_ok=True)

    corrfile1 = pickle.load(lz4.frame.open(correctionname1))
    corrfile2 = pickle.load(lz4.frame.open(correctionname2))
    proc = "Z" if "CorrZ" in correctionname1 else "W"
    var = "qT"

    # numh1 = corrfile1[proc][corrname1 + "_hist"][{"vars" : 0}]
    # numh2 = corrfile2[proc][corrname2 + "_hist"][{"vars" : 0}]
    # numh2 = corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : 0}]
    # print(corrfile1[proc])

    # print(corrfile1[proc][corrname1 + "_minnlo_ratio"][{"vars" : 0}])
    # ratio_to_minnlo = corrfile1[proc][corrname1 + "_minnlo_ratio"][{"vars" : 0}].project(var)
    nnlojet_fo = input_tools.read_nnlojet_pty_hist("/work/submit/areimers/wmass/TheoryCorrections/NNLOjet/Z/ZjNNLO/final/ptz")
    # nnlojet_fo_smooth_pt = hh.smooth_hist(nnlojet_fo.project("qT", "vars"), "qT", start_bin=4)
    scetlib_n4ll = input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_N4+0LL_combined.pkl")[{"vars" : 0}].project(var)
    scetlib_pdfas_n3ll = input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_pdf_combined.pkl")
    scetlib_pdfas_n4ll = input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl")
    scetlib_n3lo_sing = (-1*input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_n3lo_sing.pkl")[{"vars" : 0}]).project(var)

    nnlojet_fo = nnlojet_fo[{"vars" : 0}].project(var)
    # nnlojet_fo_smooth_pt = nnlojet_fo_smooth_pt[{"vars" : 0}].project(var)
    
    centers_nnlojet_fo = 0.5 * (nnlojet_fo.axes[0].edges[:-1] + nnlojet_fo.axes[0].edges[1:])
    values_nnlojet_fo = nnlojet_fo.values()
    errors_nnlojet_fo = np.sqrt(nnlojet_fo.variances())
    fit_mask = (centers_nnlojet_fo > 10) & (centers_nnlojet_fo < 60)  

    popt, pcov = curve_fit(exp_func, centers_nnlojet_fo[fit_mask], values_nnlojet_fo[fit_mask], sigma=errors_nnlojet_fo[fit_mask], absolute_sigma=True)
    A_fit, lambda_fit = popt


    # plt.figure()
    # # ratio_to_minnlo.plot()
    # mplhep.histplot(numh1, yerr=False)
    # mplhep.histplot(numh2, yerr=False)
    # plt.xlabel("qT [GeV]")
    # plt.ylabel("Correction")
    # plt.grid(True)
    # plt.savefig(f"{outfoldername}/corr1_vs_corr2_numerator_{var}.pdf")
    # plt.close()

    print(input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl"))
    fig = plot_tools.makePlotWithRatioToRef(
        hists=[
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : 0}].project(var),
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),
            input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl")[{"vars" :  "pdfCT18ZNNLO_as_0118"}].project(var),
            input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl")[{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl")[{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),

        ],
        labels=[
             "from corr: as_0118",
             "from corr: as_0116",
             "from corr: as_0120",
             "from input: as_0118",
             "from input: as_0116",
             "from input: as_0120",
            ],
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:6],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="SCETlib prediction",
        rlabel="x/0.118",
        rrange=[0.9, 1.1],
        nlegcols=1,
        xlim=None, binwnorm=1.0, baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    fig.savefig(f"{outfoldername}/corr1_vs_corr2_numerator_{var}.pdf")

    fig = plot_tools.makePlotWithRatioToRef(
        hists=[
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : 0}].project("absY"),
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : "pdfCT18ZNNLO_as_0116"}].project("absY"),
            corrfile2[proc]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ProperPdfAs_pdfas_hist"][{"vars" : "pdfCT18ZNNLO_as_0120"}].project("absY"),
        ],
        labels=[
             "pdfCT18ZNNLO_as_0118",
             "pdfCT18ZNNLO_as_0116",
             "pdfCT18ZNNLO_as_0120",
            ],
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:3],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="SCETlib prediction",
        rlabel="x/0.118",
        rrange=[0.9, 1.1],
        nlegcols=1,
        xlim=None, binwnorm=1.0, baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    fig.savefig(f"{outfoldername}/corr1_vs_corr2_numerator_absY.pdf")





    # corr_hist = corrfile1[proc][corrname1 + "_hist"][{"vars" : 0}].project(var)

    # print(nnlojet_fo)
    plt.figure()
    # mplhep.histplot(corr_hist, label="N3LO + N4LL", yerr=True)
    mplhep.histplot(nnlojet_fo, label="N3LO FO", yerr=True)
    # mplhep.histplot(nnlojet_fo_smooth_pt, label="N3LO FO (smoothed)", yerr=True)
    mplhep.histplot(scetlib_n4ll, label="N4LL resumm", yerr=True)
    mplhep.histplot(scetlib_n3lo_sing, label="N3LO sing", yerr=True)
    # mplhep.histplot(scetlib_pdfas_n3ll, label="N3LL pdfas", yerr=True)
    # mplhep.histplot(scetlib_pdfas_n4ll, label="N4LL pdfas", yerr=True)
    # xfit = np.linspace(centers_nnlojet_fo[0], centers_nnlojet_fo[-1], 500)
    # plt.plot(xfit, exp_func(xfit, *popt), label=f"Fit: A exp(-x / {lambda_fit:.1f})", color="red")
    plt.xlabel("qT [GeV]")
    plt.ylabel(f"{corrname1}")
    plt.yscale("log")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{outfoldername}/{plotoutname_base}_hist_{var}.pdf")
    plt.close()


    print(input_tools.read_scetlib_hist("/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl"))

    fig = plot_tools.makePlotWithRatioToRef(
        hists=[
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0118"}].project(var),
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0118"}].project(var),
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),
        ],
        labels=[
             "N3LL as 0.118",
             "N3LL as 0.116",
             "N3LL as 0.120",
             "N4LL as 0.118",
             "N4LL as 0.116",
             "N4LL as 0.120",
            ],
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:6],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="SCETlib prediction",
        rlabel="x/N3LL as 0.118",
        rrange=[0.9, 1.1],
        nlegcols=1,
        xlim=None, binwnorm=1.0, baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    fig.savefig(f"{outfoldername}/N3LL_N4LL_pdfas_comparison_{var}.pdf")

    fig = plot_tools.makePlotWithRatioToRef(
        hists=[
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0118"}].project(var),
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            scetlib_pdfas_n3ll[{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),
        ],
        labels=[
             "N3LL as 0.118",
             "N3LL as 0.116",
             "N3LL as 0.120",
            ],
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:3],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="SCETlib prediction",
        rlabel="x/as 0.118",
        rrange=[0.9, 1.1],
        nlegcols=1,
        xlim=None, binwnorm=1.0, baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    fig.savefig(f"{outfoldername}/N3LL_pdfas_{var}.pdf")

    fig = plot_tools.makePlotWithRatioToRef(
        hists=[
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0118"}].project(var),
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0116"}].project(var),
            scetlib_pdfas_n4ll[{"vars" : "pdfCT18ZNNLO_as_0120"}].project(var),
        ],
        labels=[
             "N4LL as 0.118",
             "N4LL as 0.116",
             "N4LL as 0.120",
            ],
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:3],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="SCETlib prediction",
        rlabel="x/as 0.118",
        rrange=[0.9, 1.1],
        nlegcols=1,
        xlim=None, binwnorm=1.0, baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    fig.savefig(f"{outfoldername}/N4LL_pdfas_{var}.pdf")


    # plt.figure()
    # mplhep.histplot(scetlib_pdfas_n3ll, label="N3LL pdfas", yerr=True)
    # mplhep.histplot(scetlib_pdfas_n4ll, label="N4LL pdfas", yerr=True)
    # plt.xlabel("qT [GeV]")
    # plt.ylabel(f"SCETlib prediction")
    # # plt.yscale("log")
    # plt.grid(True)
    # plt.legend()
    # plt.savefig(f"{outfoldername}/N3LL_N4LL_pdfas_comparison_{var}.pdf")
    # plt.close()

    # minnlo_hist = corrfile1[proc]["minnlo_ref_hist"]
    # plt.figure()
    # minnlo_hist.project(var).plot()
    # plt.xlabel("qT [GeV]")
    # plt.ylabel(f"MiNNLO")
    # plt.grid(True)
    # plt.savefig(f"{outfoldername}/{plotoutname_base}_MiNNLO_hist_{var}.pdf")
    # plt.close()

    


    # fig = plot_tools.makePlotWithRatioToRef(
    #     hists=[
    #         nnlojet_fo_smooth_pt,
    #         nnlojet_fo,
    #     ],
    #     labels=[
    #          "N3LO FO Smoothed pT only",
    #          "N3LO FO Unsmoothed",
    #         ],
    #     colors=[
    #          "black",
    #          "purple",
    #         ],
    #     xlabel="q$_{T}$ (GeV)", 
    #     ylabel="Events/bin",
    #     rlabel="x/smoothed",
    #     rrange=[0.9, 1.1],
    #     nlegcols=1,
    #     xlim=None, binwnorm=1.0, baseline=True, 
    #     yerr=True,
    #     yerr_ratio=True,
    #     linewidth=1,
    # )
    # fig.savefig(f"{outfoldername}/{plotoutname_base}_N3LO_FO_{var}_smothings.pdf")
    
    

if __name__ == "__main__":
    main()