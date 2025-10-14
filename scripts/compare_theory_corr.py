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
        "scetlib_dyturbo_ResumSelfCorrZ.pkl.lz4", 

    ]

    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    plotoutname_base = f"corrections_comparison"


    corrnames_translation = {
        "scetlib_dyturboCorrZ.pkl.lz4": "N$^{3{+}0}$LL+NNLO (nominal)",
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
    }

    # os.makedirs(args.outfolder, exist_ok=True)

    corrfiles = [pickle.load(lz4.frame.open(os.path.join(corrections_folder, c))) for c in corrections]
    proc = "Z"
    corrnames = [corr_name(c) for c in corrections]
    corrnames_pretty = [corrnames_translation[c] for c in corrections]
    var = "qT"


    for cname, cfile in zip(corrnames, corrfiles):
        print(f"\n\n--> Now at correction {cname}")
        # print(cfile[proc])
        h = cfile[proc][cname+"_hist"][{"vars": 0}]
        print(h)
        relstatunc = math.sqrt(h.sum().variance) / h.sum().value
        print(f"Relative stat unc of correction: {relstatunc}")
        h_minnlo = cfile[proc]["minnlo_ref_hist"]
        relstatunc_minnlo = math.sqrt(h_minnlo.sum().variance) / h_minnlo.sum().value
        print(f"Relative stat unc of MiNNLO reference: {relstatunc_minnlo}")
        sf = relstatunc / relstatunc_minnlo
        print(f"Should increase the MiNNLO stat unc by a factor of {sf:.02f} to reflect the stat unc in the correction")

    

    colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)]
    fig = plot_tools.makePlotWithRatioToRef(
        hists=[cfile[proc][cname+"_hist"][{"vars": 0}].project(var) for cname, cfile in zip(corrnames, corrfiles)],
        labels=corrnames_pretty,
        colors=colors[:len(corrfiles)],
        xlabel="q$_{T}$ (GeV)", 
        ylabel="Events / bin",
        rlabel="var. / nom.",
        rrange=[0.99, 1.01],
        nlegcols=1,
        xlim=None, 
        binwnorm=1.0, 
        baseline=True, 
        yerr=True,
        yerr_ratio=True,
        ratio_legend=False,
        linewidth=1,
    )
    plotoutname = f"{outfoldername}/{plotoutname_base}.pdf"
    fig.savefig(plotoutname)
    print(f"--> Wrote plot to {plotoutname}")

        

if __name__ == "__main__":
    main()