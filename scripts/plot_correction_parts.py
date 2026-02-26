import os
import pickle
import re

import lz4.frame # type: ignore

from wums import plot_tools
import matplotlib.pyplot as plt # type: ignore
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors # type: ignore
import hist # type: ignore
from copy import deepcopy



def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections"
    
    # filename_resum = ("SCETlib/ct18z_newnps_n3+0ll_lattice_pdfvars/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_withnpvars_pdfas_pdf_combined.pkl", "N$^{3{+}0}$LL")
    # filename_fosing = ("SCETlib/ct18z_nplambda_pdfvars/inclusive_Z_CT18Z_nplambda_pdfas_nnlo_sing_pdf_combined.pkl", "NNLO FO sing.")
    # filename_fo = ("DYTURBO/nnlo-scetlibmatch/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt", "NNLO FO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboN3p0LL_LatticeNP_pdfasCorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    
    # filename_resum = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_combined.pkl", "N$^{3{+}0}$LL")
    # filename_fosing = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine_nnlo_sing/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_nnlo_sing_combined.pkl", "NNLO FO sing.")
    # filename_fo = ("DYTURBO/nnlo-scetlibmatch-13TeV-CT18Z-finer-bin/scalevariations/z0/results_z-2d-nnlo-vj-CT18ZNNLO-{scale}-scetlibmatch.txt", "NNLO FO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    # # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo__LatticeNP_CT18Z_N3p0LL_N2L0_test_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    
    # filename_resum = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine_pdfvars/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_pdfas_pdf_combined.pkl", "N$^{3{+}0}$LL")
    # filename_fosing = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine_pdfvars_nnlo_sing/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_pdfas_nnlo_sing_pdf_combined.pkl", "NNLO FO sing.")
    # filename_fo = ("DYTURBO/nnlo-scetlibmatch-13TeV-CT18Z-finer-bin/pdfvariations/CT18ZNNLO_as/z0/results_z-2d-nnlo-vj-CT18ZNNLO_as-member{i}-scetlibmatch.txt", "NNLO FO")
    # # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_pdfas_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo__LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfas_test_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    
    # filename_resum = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine_pdfvars/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_pdfvars_combined.pkl", "N$^{3{+}0}$LL")
    # filename_fosing = ("SCETlib/com13_ct18z_newnps_n3+0ll_lattice_fine_pdfvars_nnlo_sing/inclusive_Z_COM13_CT18Z_N3+0LL_lattice_fine_pdfvars_nnlo_sing_combined.pkl", "NNLO FO sing.")
    # filename_fo = ("DYTURBO/nnlo-scetlibmatch-13TeV-CT18Z-finer-bin/pdfvariations/CT18ZNNLO/z0/results_z-2d-nnlo-vj-CT18ZNNLO-member{i}-scetlibmatch.txt", "NNLO FO")
    # # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_pdfvars_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo__LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfvars_test_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    
    # filename_resum = ("SCETlib/com13_nnpdf31_newnps_n3+0ll_lattice_fine_pdfvars/inclusive_Z_COM13_NNPDF31_N3+0LL_lattice_fine_pdfas_combined.pkl", "N$^{3{+}0}$LL")
    # filename_fosing = ("SCETlib/com13_nnpdf31_newnps_n3+0ll_lattice_fine_pdfvars_nnlo_sing/inclusive_Z_COM13_NNPDF31_N3+0LL_lattice_fine_pdfas_nnlo_sing_combined.pkl", "NNLO FO sing.")
    # filename_fo = ("DYTURBO/nnlo-scetlibmatch-13TeV-NNPDF31-finer-bin/pdfvariations/NNPDF31_nnlo_hessian_pdfas_as/z0/results_z-2d-nnlo-vj-NNPDF31_as-member{i}-scetlibmatch.txt", "NNLO FO")
    # filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_NNPDF31_N3p0LL_N2L0_pdfas_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")
    
    filename_resum = ("SCETlib/com13_pdf4lhc21_newnps_n3+0ll_lattice_fine_pdfvars/inclusive_Z_COM13_PDF4LHC21_N3+0LL_lattice_fine_pdfas_combined.pkl", "N$^{3{+}0}$LL")
    filename_fosing = ("SCETlib/com13_pdf4lhc21_newnps_n3+0ll_lattice_fine_pdfvars_nnlo_sing/inclusive_Z_COM13_PDF4LHC21_N3+0LL_lattice_fine_pdfas_nnlo_sing_combined.pkl", "NNLO FO sing.")
    filename_fo = ("DYTURBO/nnlo-scetlibmatch-13TeV-PDF4LHC21-finer-bin/pdfvariations/PDF4LHC21_40_pdfas_as/z0/results_z-2d-nnlo-vj-PDF4LHC21_40_pdfas_as_member{i}-mur1-muf1-scetlibmatch.txt", "NNLO FO")
    filename_fullcorr = ("/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_PDF4LHC21_N3p0LL_N2L0_pdfas_CorrZ.pkl.lz4", "Full N$^{3{+}0}$LL+NNLO")

    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    os.makedirs(outfoldername, exist_ok=True)

    # variations = [[{"vars" :  "pdfCT18ZNNLO_as_0118", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.118\\pm0.002$)"], [{"vars" :  "pdfCT18ZNNLO_as_0116", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.116$)"], [{"vars" :  "pdfCT18ZNNLO_as_0120", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.120$)"]]
    # vartag = "asvars"

    # variations = [[{"vars" :  "pdf0", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.118\\pm0.002$)"], [{"vars" :  "pdf101", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.116$)"], [{"vars" :  "pdf102", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.120$)"]]
    # vartag = "asvars"

    variations = [[{"vars" :  "pdf0", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.118\\pm0.001$)"], [{"vars" :  "pdf41", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.117$)"], [{"vars" :  "pdf42", "Q": 0}, "($\\alpha_{\\mathrm{s}} = 0.119$)"]]
    vartag = "asvars"

    # variations = [[{"vars" :  "pdf0", "Q": 0}, "(PDF vars)"], [{"vars" :  "pdf1", "Q": 0}, "(PDF 1)"], [{"vars" :  "pdf2", "Q": 0}, "(PDF 2)"], [{"vars" :  "pdf3", "Q": 0}, "(PDF 3)"], [{"vars" :  "pdf4", "Q": 0}, "(PDF 4)"]]
    # vartag = "pdfvars"

    # variations = [({"vars" :  0, "Q": 0}, "($\\mu_{r} = 1$)"), ({"vars" :  "kappaFO0.5-kappaf2.", "Q": 0}, "($\\mu_{r} = 0.5$)"), ({"vars" :  "kappaFO2.-kappaf0.5", "Q": 0}, "($\\mu_{r} = 2$)")]
    # vartag = "murvars"

    # variations = [[{"vars" :  0, "Q": 0}, "nominal"]]
    # vartag = "nominal"

    # variations = [({"vars" :  "kappaFO0.5-kappaf2.", "Q": 0}, "($\\mu_{r} = 0.5$)")]
    # vartag = "mur0p5"

    variations_parts = []
    for v in variations:
        v_part = deepcopy(v)
        v_part_dict = v_part[0]
        v_part_dict["Q"] = 1
        v_part[0] = v_part_dict
        variations_parts.append(v_part)


    plot_tag = os.path.basename(filename_fullcorr[0]).split(".")[0]
    
    hist_resum = input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum[0]))
    hist_fosing = input_tools.read_scetlib_hist(os.path.join(genpath, filename_fosing[0]))
    hist_fullcorr = pickle.load(lz4.frame.open(filename_fullcorr[0]))["Z"][f"{plot_tag.replace("CorrZ", "")}hist"]
    print(hist_resum.axes)
    print(hist_fullcorr.axes)
        
    if filename_fo[0].startswith("DYTURBO"):
        # hist_fo = input_tools.read_dyturbo_vars_hist(os.path.join(genpath, filename_fo[0]), var_axis=hist.axis.StrCategory(['pdfCT18ZNNLO_as_0118', 'pdfCT18ZNNLO_as_0116', 'pdfCT18ZNNLO_as_0120'],     name="vars"), axes=['Y', 'qT'], charge=0)
        hist_fo = input_tools.read_dyturbo_vars_hist(os.path.join(genpath, filename_fo[0]), var_axis=hist_fosing.axes["vars"], axes=["Q", 'Y', 'qT'], charge=0)
    elif filename_fo[0].startswith("NNLOjet"):
        y_axes = np.array([-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ])
        hist_fo =input_tools.read_nnlojet_pty_hist(os.path.join(genpath, filename_fo[0]), ybins=y_axes, charge=0)

    
    vars_resum      = [hist_resum[v[0]].project("qT") for v in variations_parts]
    vars_fo         = [hist_fo[v[0]].project("qT") for v in variations_parts]
    vars_fullcorr   = [hist_fullcorr[v[0]].project("qT") for v in variations]
    vars_fosing_neg = [hist_fosing[v[0]].copy() for v in variations_parts]
    for v in vars_fosing_neg:
        v.view(flow=True)[...] *= -1
    vars_fosing_neg = [h.project("qT") for h in vars_fosing_neg]

    hists = vars_fullcorr+vars_resum+vars_fosing_neg+vars_fo
    labels = [f"{filename_fullcorr[1]} {variations[0][1]}"] + [""]*(len(variations)-1) + [f"{filename_resum[1]} {variations_parts[0][1]}"] + [""]*(len(variations_parts)-1) + [f"{filename_fosing[1]} {variations_parts[0][1]}"] + [""]*(len(variations_parts)-1) + [f"{filename_fo[1]} {variations_parts[0][1]}"] + [""]*(len(variations_parts)-1)
    linestyles = (["-"] + [(0, (1, 0.75))]*(len(variations)-1))*4
    colors = []
    for i in range(4):
        colors += [f"{mcolors.to_hex(list(plt.get_cmap('tab10').colors)[i])}"] * len(variations)

    fig = plot_tools.makePlotWithRatioToRef(
        hists=hists,
        labels=labels,
        colors=colors,
        linestyles=linestyles,
        xlabel="$q_{T}$ (GeV)", 
        ylabel="Prediction (pb/GeV)",
        # ylabel="Prediction (pb/bin)",
        rlabel="x / full pred.",
        rrange=[0.8, 1.2],
        nlegcols=1,
        binwnorm=1.0, 
        baseline=True, 
        yerr=False,
        ratio_legend=False,
        linewidth=1.3,
        logy=False,
        legtext_size=20,
        xlim=[0., 100.],
    )
    # fig.text(0.2, 0.88, "$Q \\in [10, 60]$", ha="left", va="bottom", fontsize=20)
    fig.text(0.2, 0.88, "$Q \\in [60, 120]$", ha="left", va="bottom", fontsize=20)
    outfilename = f"{outfoldername}/parts_{plot_tag}_{vartag}_qT.pdf"
    fig.savefig(outfilename)
    print(f"--> Saved plot to {outfilename}")

    # for variation in variations:
    #     var_resum      = hist_resum[variation]
    #     var_fosing     = hist_fosing[variation]
    #     var_fosing_neg = var_fosing.copy()  # make a copy
    #     var_fosing_neg.view(flow=True)[...] *= -1
    #     var_fo         = hist_fo[variation]
    #     var_fullcorr   = hist_fullcorr[variation]
    
    #     fig = plot_tools.makePlotWithRatioToRef(
    #         hists=[var_fullcorr.project("qT"), var_resum.project("qT"), var_fosing_neg.project("qT"), var_fo.project("qT")],
    #         labels=[filename_fullcorr[1], filename_resum[1], filename_fosing[1], filename_fo[1]],
    #         colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)][:4],
    #         xlabel="$q_{T}$ (GeV)", 
    #         ylabel="Prediction",
    #         rlabel="x/full corr.",
    #         rrange=[0.5, 1.5],
    #         nlegcols=1,
    #         xlim=None, binwnorm=1.0, 
    #         baseline=True, 
    #         yerr=True,
    #         ratio_legend=False,
    #         linewidth=1,
    #         logy=False,
    #     )
    #     outfilename = f"{outfoldername}/parts_{plot_tag}_{variation["vars"]}_qT.pdf"
    #     fig.savefig(outfilename)
    #     print(f"--> Saved plot to {outfilename}")

if __name__ == "__main__":
    main()