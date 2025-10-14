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
from itertools import product


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
    
    filenames_resum = ["SCETlib/ct18z_nplambda_pdfvars_arnetest/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl", "SCETlib/ct18z_nplambda_pdfvars_arnetest/inclusive_Z_CT18Z_nplambda_N3p1LL_pdfas_pdf_combined.pkl"]
    # filenames_resum = ["SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_pdf_combined.pkl"]
    
    filenames_fosing = ["SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl", "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl", "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl"]

    filenames_fullcorr = ["/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfasCorrZ.pkl.lz4", "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_ResumSelf_pdfasCorrZ.pkl.lz4"]
    # filenames_fullcorr = ["/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_pdfasCorrZ.pkl.lz4"]
    # filenames_fullcorr = ["/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfasCorrZ.pkl.lz4"]

    filenames_fo = ["DYTURBO/nnlo-scetlibmatch/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt", "DYTURBO/nnlo-scetlibmatch/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt", "DYTURBO/nnlo-scetlibmatch/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt"]


    vars_to_plot = [{"vars" :  "pdfCT18ZNNLO_as_0118"}, {"vars" :  "pdfCT18ZNNLO_as_0116"}, {"vars" :  "pdfCT18ZNNLO_as_0120"}]
    # vars_to_plot = [{"vars" :  0}, {"vars" :  1}, {"vars" :  2}]
    plot_tag_resum = "_asFromN4LL"
    # plot_tag_resum = "_asFromN3LL"
    plot_tag_fosing = "_asFromN2LO"
    plot_tag_fullcorr = "_asFromN4LL"
    # plot_tag_fullcorr = "_asFromN3LL"
    # plot_tag_fullcorr = "N4LLN2LO_asFromN4LL"
    plot_tag_fo = "_asFromN2LO"

    hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, fnr)) for fnr in filenames_resum]
    hists_fosing = [input_tools.read_scetlib_hist(os.path.join(genpath, fns)) for fns in filenames_fosing]
    hists_fo = [input_tools.read_dyturbo_vars_hist(os.path.join(genpath, fnf), var_axis=hist.axis.StrCategory(['pdfCT18ZNNLO_as_0118', 'pdfCT18ZNNLO_as_0116', 'pdfCT18ZNNLO_as_0120'], name="vars"), axes=['Y', 'qT'], charge=0) for fnf in filenames_fo]

    histnames_fullcorr = [f"{os.path.basename(fnf).replace('CorrZ.pkl.lz4', '')}_hist" for fnf in filenames_fullcorr]
    histnames_fullcorr = [f"{hnf.replace('N3LO', 'N2LO').replace('nnlojet', 'dyturbo')}" for hnf in histnames_fullcorr]
    hists_fullcorr = [pickle.load(lz4.frame.open(fnf))["Z"][hnf] for (fnf, hnf) in zip(filenames_fullcorr, histnames_fullcorr)]
    # hists_fullcorr = [pickle.load(lz4.frame.open(fnf))["Z"]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfas_hist"] for fnf in filenames_fullcorr]
    # hists_fullcorr = [pickle.load(lz4.frame.open(fnf))["Z"]["scetlib_dyturboCT18Z_pdfas_hist"] for fnf in filenames_fullcorr]

    vars_resum  = [h[v] for (h, v) in product(hists_resum, vars_to_plot)]
    vars_fosing  = [h[v] for (h, v) in zip(hists_fosing, vars_to_plot)]
    vars_fullcorr  = [h[v] for (h, v) in product(hists_fullcorr, vars_to_plot)]
    vars_fo  = [h[v] for (h, v) in zip(hists_fo, vars_to_plot)]

    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    os.makedirs(outfoldername, exist_ok=True)


    variables = ["qT", "Y"]

    for v in variables:
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[var_resum.project(v) for var_resum in vars_resum],
            labels=[str(var_to_plot["vars"]) for fnr in filenames_resum for var_to_plot in vars_to_plot],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(vars_resum)],
            xlabel=v, 
            ylabel="Resumm. prediction",
            rlabel="ratio",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_ascomparison_resumm{plot_tag_resum}_{v}.pdf")

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[var_fosing.project(v) for var_fosing in vars_fosing],
            labels=[str(var_to_plot["vars"]) for var_to_plot in vars_to_plot],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(vars_to_plot)],
            xlabel=v, 
            ylabel="FO sing. prediction",
            rlabel="ratio",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_ascomparison_FOsingular{plot_tag_fosing}_{v}.pdf")

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[var_fo.project(v) for var_fo in vars_fo],
            labels=[str(var_to_plot["vars"]) for var_to_plot in vars_to_plot],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(vars_to_plot)],
            xlabel=v, 
            ylabel="FO. prediction",
            rlabel="ratio",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_ascomparison_FO{plot_tag_fo}_{v}.pdf")

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[var_fullcorr.project("absY" if v == "Y" else v) for var_fullcorr in vars_fullcorr],
            labels=[str(var_to_plot["vars"]) for fnf in filenames_fullcorr for var_to_plot in vars_to_plot],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:len(vars_fullcorr)],
            xlabel=v, 
            ylabel="Full combined prediction",
            rlabel="ratio",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            legtext_size=15,
            logy=True,
        )
        fig.savefig(f"{outfoldername}/correction_ascomparison_FullCorr{plot_tag_fullcorr}_{v}.pdf")

if __name__ == "__main__":
    main()