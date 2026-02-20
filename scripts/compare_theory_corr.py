import argparse
import os
import pickle
import re
import math

import lz4.frame

from utilities import common
from wremnants import theory_corrections
from wums import boostHistHelpers as hh
from wums import logging, output_tools, plot_tools
from scripts.corrections.make_theory_corr import read_corr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from utilities.io_tools import input_tools
import mplhep
import numpy as np
from scipy.optimize import curve_fit
from hist import tag


def corr_name(corrf):
    if not corrf.endswith(".pkl.lz4"):
        raise ValueError(f"File {corrf} is not a lz4 compressed pickle file")

    match = re.match(r"(.*)Corr[Z|W]\.pkl\.lz4", corrf)
    return match[1]



def main():
    
    corrections_folder = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections"
    # corrections = ["scetlib_dyturboCorrZ.pkl.lz4", "scetlib_dyturboN3p1LLCorrZ.pkl.lz4", "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4", "scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4", "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4"]
    corrections = [
        "scetlib_dyturboCorrZ.pkl.lz4", 
        # "scetlib_dyturboN3p1LLCorrZ.pkl.lz4", 
        # "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4", 
        # "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixed_N2p1LLFixedCorrZ.pkl.lz4",
        # "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixed_N2p1LLFixed_ResumSelfCorrZ.pkl.lz4",
        # "scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4",
        # "scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ResumSelfCorrZ.pkl.lz4",
        # "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4",
        # "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_ResumSelfCorrZ.pkl.lz4",
        # "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_AsAtN4LL_ResumSelfCorrZ.pkl.lz4",
        # "scetlib_dyturbo_ResumSelfCorrZ.pkl.lz4", 
        # "scetlib_dyturbo_NewNPModel_LatticeValsAndVarsCorrZ.pkl.lz4",
        # "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_WithLatticeFOSingCorrZ.pkl.lz4",
        # "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun_pdfasCorrZ.pkl.lz4",
        "scetlib_dyturboN3p0LL_LatticeNPCorrZ.pkl.lz4",
        # "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorrZ.pkl.lz4",
        # "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_CorrZ.pkl.lz4",

    ]

    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    plotoutname_base = f"corrections_comparison"


    corrnames_translation = {
        "scetlib_dyturboCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (old NP model)",
        "scetlib_dyturboN3p1LLCorrZ.pkl.lz4": "N$^{2{+}1}$LL+NNLO",
        "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4": "N$^{2{+}1}$LL+N$^{3}$LO",
        "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixed_N2p1LLFixedCorrZ.pkl.lz4": "N$^{3{+}1}$LL+N$^{3}$LO",
        "scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixed_N2p1LLFixed_ResumSelfCorrZ.pkl.lz4": "N$^{3{+}1}$LL+N$^{3}$LO (resum. from Arne)",
        "scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4": "N$^{4{+}0}$LL+NNLO",
        "scetlib_dyturbo_N4p0LLN2LOUnsmoothed_N3pLLFixed_ResumSelfCorrZ.pkl.lz4": "N$^{4{+}0}$LL+NNLO (resum. from Arne)",
        "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorrZ.pkl.lz4": "N$^{4{+}0}$LL+N$^{3}$LO",
        "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_ResumSelfCorrZ.pkl.lz4": "N$^{4{+}0}$LL+N$^{3}$LO (resum. from Arne)",
        "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixed_AsAtN4LL_ResumSelfCorrZ.pkl.lz4": "N$^{4{+}0}$LL+N$^{3}$LO (as@N$^{4}$4LL, resum. from Arne)",
        "scetlib_dyturbo_ResumSelfCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (remade N$^{3{+}0}$LL piece)",
        "scetlib_dyturbo_NewNPModel_LatticeValsAndVarsCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars)",
        "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_WithLatticeFOSingCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (lattice vals/vars, +lattice FO singular)",
        "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun_pdfasCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (lattice vals/vars, rerun pdfas)",
        "scetlib_dyturboN3p0LL_LatticeNPCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (new NP model)",
        "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (lattice, pdfas file)",
        "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_CorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (new NP, finer binning)",
    }

    vars_translation = {
        "qT": "q$_{T}$ (GeV)",
        "absY": "|Y|",
    }

    yaxis_per_var = {
        "qT": "Prediction (pb / GeV)",
        "absY": "Prediction (pb / bin)",
    }

    ratio_ranges_per_var = {
        "qT": [0.90, 1.10],
        "absY": [0.995, 1.005],
    }

    # os.makedirs(args.outfolder, exist_ok=True)

    corrfiles = [pickle.load(lz4.frame.open(os.path.join(corrections_folder, c))) for c in corrections]
    proc = "Z"
    corrnames = [corr_name(c) for c in corrections]
    corrnames_pretty = [corrnames_translation[c] for c in corrections]
    vars = ["qT", "absY"]

    for var in vars:
        colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)]
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[cfile[proc][cname+"_hist"][{"vars": 0}].project(var) for cname, cfile in zip(corrnames, corrfiles)],
            labels=corrnames_pretty,
            colors=colors[:len(corrfiles)],
            xlabel=vars_translation[var], 
            ylabel=yaxis_per_var[var],
            rlabel="new / old",
            rrange=ratio_ranges_per_var[var],
            nlegcols=1,
            xlim=None, 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            # yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
        )
        plotoutname = f"{outfoldername}/{plotoutname_base}_{var}.pdf"
        fig.savefig(plotoutname)
        print(f"--> Wrote plot to {plotoutname}")

        

if __name__ == "__main__":
    main()