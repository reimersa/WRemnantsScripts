import os
import pickle

import lz4.frame

from wums import boostHistHelpers as hh
from wums import output_tools, plot_tools
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import mplhep
import numpy as np
from scipy.optimize import curve_fit

import h5py




def main():
    plot_variations = True

    # corrname = "scetlib_nnlojet_N3p0LLN2LOUnsmoothedCorr"
    # corrname = "scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorr"
    # corrname = "scetlib_nnlojetN4p0LLN3LOUnsmoothedCorr"
    # corrname = "scetlib_dyturboCorr"
    # corrname = "scetlib_dyturbo_NewNPModel_OldValsAndVarsCorr"
    # corrname = "scetlib_dyturbo_NewNPModel_LatticeValsAndVarsCorr"
    corrname = "scetlib_dyturboN3p0LL_LatticeNPCorr"
    # corrname = ""

    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/mz_dilepton_{corrname}_maxFiles_m1.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/mz_dilepton_{corrname}_maxFiles_1__perBinStatCorrUnc.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/mz_dilepton_{corrname}_maxFiles_1_perBinStatCorrUncNdim.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/mw_with_mu_eta_pt_{corrname}_maxFiles_m1_FixedAlphaS.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/mz_dilepton_{corrname}_maxFiles_m1_FixedAlphaS.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_{corrname}_maxFiles_m1_ct18z_pdfByHelicity.hdf5"
    # infilename = f"/work/submit/areimers/wmass/WRemnants/wremnants-data/data/PDFs/w_z_gen_dists_maxFiles_m1_ct18z_pdfByHelicity_skimmed.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/mz_dilepton_{corrname}_maxFiles_m1_FixedAlphaS_AsAndPdfFromHels_Inputs4d_OwnPDFVars.hdf5"
    infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/mz_dilepton_{corrname}_maxFiles_m1_FixedAlphaS_AsAndPdfFromHels_Inputs4d_OwnPDFVars_InclOldNPPDFAsVars.hdf5"
    
    

    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/variations/{corrname}"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/variations/{corrname}_perBinStatCorrUnc"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/variations/{corrname}_perBinStatCorrUncNdim"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/W/variations/{corrname}_FixedAlphaS"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/ForAlphaS/WRemDev/variations/{corrname}_FixedAlphaS"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/ForAlphaS/WRemDev/variations/{corrname}_FixedAlphaS_AsAndPdfFromHels_Inputs4d_OwnPDFVars"
    outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/ForAlphaS/WRemDev/variations/{corrname}_FixedAlphaS_AsAndPdfFromHels_Inputs4d_OwnPDFVars_InclOldNPPDFAsVars"
    

    with h5py.File(infilename, "r") as h5file:
        results = input_tools.load_results_h5py(h5file)

    print(results.keys())

    samplename = "ZmumuPostVFP"
    # samplename = "WminusmunuPostVFP"

    print(results[samplename]["output"].keys())



    # histname = f"nominal_ptll_{corrname}"
    # histname = f"nominal_{corrname}"
    # histname = f"nominal_{corrname.replace("Corr", "FlavDepNP")}"
    # histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr"
    # histname = "nominal_gen_pdfCT18Z"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr"
    histname = "nominal_scetlib_dyturboCT18ZVarsCorr"
    # histname = "nominal_pdfAlphaSByHelicity"
    hist = input_tools.read_and_scale(infilename, samplename, histname)
    print(hist)
    
    # raise ValueError()

    varaxisname = "vars"
    project_on = "ptll"
    # vars = ["pdf14", "pdf26", "pdf27"]
    vars = ["pdf14", "pdf26", "pdf27"]
    # vars = list(hist.axes)
    # print(vars)



    if plot_variations:
        os.makedirs(outfolder, exist_ok=True)
        # nominal = hist[{"vars" : 0}]
        nominal = hist[{varaxisname : 0}].project(project_on)
        print(nominal)
        print(nominal.axes)
        # raise ValueError("stop")
        # for var in vars:
        #     varhist = hist[{varaxisname : var}]

        #     fig = plot_tools.makePlotWithRatioToRef(
        #         hists=[
        #             nominal,
        #             varhist,
        #         ],
        #         labels=[
        #              "Nominal",
        #              f"{var}",
        #             ],
        #         colors=[
        #              "black",
        #              "red",
        #             ],
        #         xlabel="p$_{T}^{\ell\ell}$ (GeV)", 
        #         ylabel="Events/bin",
        #         rlabel="x/nominal",
        #         rrange=[0.9, 1.1],
        #         nlegcols=1,
        #         xlim=None, binwnorm=1.0, baseline=True, 
        #         ratio_legend=False,
        #         yerr=False,
        #         yerr_ratio=False,
        #         linewidth=1,
        #     )
        #     fig.savefig(f"{outfolder}/{var}.pdf")
        #     plt.close(fig)

        colors = ["red", "blue", "green", "orange", "pink", "violet"]
        for idxplot in range(0, len(vars), 6):
            varnames = vars[idxplot:idxplot+6]
            varhists = [hist[{varaxisname : v}].project(project_on) for v in varnames]
            fig = plot_tools.makePlotWithRatioToRef(
                hists=[nominal] + varhists,
                labels=["Nominal"] + varnames,
                colors=["black"] + colors[0:len(varnames)],
                xlabel="p$_{T}^{V}$ (GeV)", 
                ylabel="Events/bin",
                rlabel="x/nominal",
                rrange=[0.95, 1.05],
                nlegcols=1,
                # xlim=None, 
                xlim=[0, 130.], 
                binwnorm=1.0, 
                baseline=True, 
                ratio_legend=False,
                yerr=False,
                # yerr_ratio=False,
                linewidth=1,
            )
            outfilename = f"{outfolder}/summary_{int(idxplot/5):02d}.pdf"
            fig.savefig(outfilename)
            print(f"--> Saved {outfilename}")
            plt.close(fig)

    
    

if __name__ == "__main__":
    main()