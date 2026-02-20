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

    # corrname = "scetlib_dyturboCorr"
    # corrname = "scetlib_dyturboN3p0LL_LatticeNPCorr"
    # corrname = "scetlib_dyturboN3p0LL_LatticeNP_RelToMinnloCorr"
    # corrname = "scetlib_dyturboMSHT20Corr"
    corrname = "scetlib_dyturboN3p0LL_LatticeNP_MSHT20Corr"

    # tag = "FixedAlphaS_OwnPDFVars"
    # tag = "FixedAlphaS_HackedPDFVarsFromHels_OwnPDFVars"
    # tag = "FixedAlphaS_AsAndPdfFromHels_Inputs4d"
    # tag = "msht20_OldNP_AllMSHT20"
    tag = "msht20_NewNP_AllMSHT20"

    # infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/mz_dilepton_{corrname}_maxFiles_m1_{tag}.hdf5"
    # outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/ForAlphaS/WRemDev/variations/{corrname}_{tag}"
    infilename = f"/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/NewSteer/mz_dilepton_{corrname}_maxFiles_m1_{tag}.hdf5"
    outfolder = f"/work/submit/areimers/wmass/plots/histmaker/Z/ForAlphaS/WRemDev/NewSteer/variations/{corrname}_{tag}"

    with h5py.File(infilename, "r") as h5file:
        results = input_tools.load_results_h5py(h5file)

    print(results.keys())

    samplename = "ZmumuPostVFP"

    print(results[samplename]["output"].keys())



    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr"
    # histname = "nominal_pdfCT18Z_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorrUncertByHelicity"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_RelToMinnlo_CT18ZVarsCorr"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorrUncertByHelicity"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_pdfasCorrByHelicity"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNPCorr"
    # histname = "nominal_scetlib_dyturboCT18Z_pdfasCorr"
    # histname = "nominal_pdfCT18ZUncertByHelicity"
    # histname = "nominal_scetlib_dyturboCT18ZVarsCorrUncertByHelicity"
    # histname = "nominal_scetlib_dyturboMSHT20_pdfasCorrByHelicity"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorrByHelicity"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorrUncertByHelicity"
    # histname = "nominal_scetlib_dyturboMSHT20VarsCorrUncertByHelicity"
    histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorrUncertByHelicity"

    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_RelToMinnlo_pdfasCorr"
    # histname = "nominal_scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr"
    
    hist = input_tools.read_and_scale(infilename, samplename, histname)
    print(hist)
    

    varaxisname = "vars"
    # varaxisname = "pdfVar"
    project_on = "ptll"

    # vars = [f"pdf0", "pdf2", "pdf5"]
    # vars = [f"pdf{idx}" for idx in range(1, 65)]
    vars = [idx for idx in range(1, 65)]
    # vars = [f"pdf{idx}" for idx in range(1, 59)]
    # vars = []
    # for idx in range(1, 30):
    #     vars += [f"pdf{idx}CT18ZDown", f"pdf{idx}CT18ZUp"]

    # vars = [f"pdfCT18ZNNLO_as_0118", "pdfCT18ZNNLO_as_0116", "pdfCT18ZNNLO_as_0120"]
    # vars = [0, 1, 2]



    if plot_variations:
        os.makedirs(outfolder, exist_ok=True)
        nominal = hist[{varaxisname : 0}].project(project_on)
        print(nominal)
        print(nominal.axes)

        colors = ["red", "blue", "green", "orange", "pink", "violet"]
        for idxplot in range(0, len(vars), 5):
            varnames = vars[idxplot:idxplot+5]
            varhists = [hist[{varaxisname : v}].project(project_on) for v in varnames]
            fig = plot_tools.makePlotWithRatioToRef(
                hists=[nominal] + varhists,
                labels=["Nominal"] + varnames,
                colors=["black"] + colors[0:len(varnames)],
                xlabel="p$_{T}^{ll}$ (GeV)", 
                ylabel="Events/bin",
                rlabel="x/nominal",
                rrange=[0.95, 1.05],
                nlegcols=1,
                xlim=None, 
                ylim=[0, 30.], 
                binwnorm=1.0, 
                baseline=True, 
                ratio_legend=False,
                yerr=False,
                # yerr_ratio=False,
                linewidth=1,
            )
            fig.savefig(f"{outfolder}/summary_{int(idxplot/5):02d}.pdf")
            plt.close(fig)

    
    

if __name__ == "__main__":
    main()