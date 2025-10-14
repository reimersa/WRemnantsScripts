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
    
    # filename_resum = "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_N4LL_pdfas_pdf_combined.pkl"
    # filename_resum = "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_pdf_combined.pkl"
    # filename_resum = ("SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_N4+0LL_combined.pkl", "N$^{4+0}$LL")
    # filename_resum = ("SCETlib/ct18z_nplambda_n3+0ll/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl", "N$^{3+0}$LL")
    filename_resum = ("SCETlib/ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_combined.pkl", "N$^{3+0}$LL (lattice)")

    # filename_fosing = "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl"
    # filename_fosing = "SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl"
    # filename_fosing = ("SCETlib/ct18z_nplambda_n4+0ll/inclusive_Z_CT18Z_nplambda_n3lo_sing.pkl", "N$^{3}$LO FO sing.")
    filename_fosing = ("SCETlib/ct18z_nplambda_scalevars/inclusive_Z_CT18Z_nplambda_scalevars_nnlo_sing_combined.pkl", "NNLO FO sing.")

    filename_fo = ("DYTURBO/nnlo-scetlibmatch/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt", "NNLO FO")
    # filename_fo = ("NNLOjet/Z/ZjNNLO/final/ptz", "N$^{3}$LO FO")
    
    # filename_fullcorr = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfasCorrZ.pkl.lz4"
    # filename_fullcorr = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_pdfasCorrZ.pkl.lz4"
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4", "N$^{4{+}0}$LL+N$^{3}$LO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboCorrZ.pkl.lz4", "N$^{3{+}0}$LL+NNLO")
    filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_NewNPModel_LatticeValsAndVarsCorrZ.pkl.lz4", "N$^{3{+}0}$LL+NNLO (lattice)")

    # variations = [{"vars" :  "pdfCT18ZNNLO_as_0118"}, {"vars" :  "pdfCT18ZNNLO_as_0116"}, {"vars" :  "pdfCT18ZNNLO_as_0120"}]
    variations = [{"vars" :  0}]
    # variations = [{"vars" :  0}, {"vars" :  1}, {"vars" :  2}]

    for variation in variations:

        # plot_tag = "_scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfasCorrZ"
        # plot_tag = "_scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_pdfasCorrZ"
        # plot_tag = "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorrZ"
        plot_tag = os.path.basename(filename_fullcorr[0]).split(".")[0]
    
        hist_resum = input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum[0]))
        hist_fosing = input_tools.read_scetlib_hist(os.path.join(genpath, filename_fosing[0]))
        
        if filename_fo[0].startswith("DYTURBO"):
            hist_fo = input_tools.read_dyturbo_vars_hist(os.path.join(genpath, filename_fo[0]), var_axis=hist.axis.StrCategory(['pdfCT18ZNNLO_as_0118', 'pdfCT18ZNNLO_as_0116', 'pdfCT18ZNNLO_as_0120'],     name="vars"), axes=['Y', 'qT'], charge=0)
        elif filename_fo[0].startswith("NNLOjet"):
            y_axes = np.array([-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ])
            hist_fo =input_tools.read_nnlojet_pty_hist(os.path.join(genpath, filename_fo[0]), ybins=y_axes, charge=0)

        # hist_fullcorr = pickle.load(lz4.frame.open(filename_fullcorr))["Z"]["scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_PdfAsFromN4LL_pdfas_hist"]
        # hist_fullcorr = pickle.load(lz4.frame.open(filename_fullcorr))["Z"]["scetlib_dyturboCT18Z_pdfas_hist"]
        hist_fullcorr = pickle.load(lz4.frame.open(filename_fullcorr[0]))["Z"][f"{plot_tag.replace("CorrZ", "")}_hist"]
    
        var_resum     = hist_resum[variation]
        var_fosing    = hist_fosing[variation]
        var_fosing_neg = var_fosing.copy()  # make a copy
        var_fosing_neg.view(flow=True)[...] *= -1
        var_fo    = hist_fo[variation]
        var_fullcorr  = hist_fullcorr[variation]
    
        outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
        os.makedirs(outfoldername, exist_ok=True)
    
    
    
            
            
        fo_nonsing = var_fo.project("qT").copy()
        fo_nonsing.view(flow=False)[...] = var_fo.project("qT").view(flow=False) - var_fosing.project("qT").view(flow=False)

        sum_of_parts = var_resum.project("qT").copy()
        sum_of_parts.view(flow=False)[...] = var_resum.project("qT").view(flow=False)[...] + var_fo.project("qT").view(flow=False)[...] + - var_fosing.project("qT").view(flow=False)

        fig = plot_tools.makePlotWithRatioToRef(
            # hists=[var_fullcorr.project("qT"), var_resum.project("qT"), var_fosing_neg.project("qT"), var_fo.project("qT"), fo_nonsing],
            # labels=["Full correction", "Resumm.", "FO sing.", "FO", "FO nonsing."],
            # colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:5],
            hists=[var_fullcorr.project("qT"), var_resum.project("qT"), var_fosing_neg.project("qT"), var_fo.project("qT")],
            labels=[filename_fullcorr[1], filename_resum[1], filename_fosing[1], filename_fo[1]],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:4],
            # hists=[var_fullcorr.project("qT"), var_resum.project("qT"), var_fosing_neg.project("qT"), var_fo.project("qT"), sum_of_parts],
            # labels=["Full correction", "Resumm.", "FO sing.", "FO", "Sum of parts"],
            # colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:5],
            xlabel="qT (GeV)", 
            ylabel="Prediction",
            rlabel="x/full corr.",
            rrange=[0.9, 1.1],
            nlegcols=1,
            xlim=None, binwnorm=1.0, baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=False,
        )
        fig.savefig(f"{outfoldername}/parts{plot_tag}_{variation["vars"]}_qT.pdf")

if __name__ == "__main__":
    main()