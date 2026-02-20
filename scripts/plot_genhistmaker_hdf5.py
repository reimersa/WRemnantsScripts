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

    # corrname = "scetlib_dyturboMSHT20_pdfasCorr"
    # corrname = "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr"
    # corrname = "scetlib_dyturboCT18Z_pdfasCorr"
    # corrname = "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr"
    # corrname = "scetlib_dyturboMSHT20VarsCorr"
    corrname = "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr"

    # infilename = f"/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrectionsByHelicity/AlphaS/w_z_gen_dists_{corrname}_maxFiles_m1_skimmed.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_scetlib_dyturboMSHT20_pdfasCorr_maxFiles_m1_msht20_msht20_alphasByHelicity_TEST.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr_maxFiles_m1_msht20_msht20_alphasByHelicity_TEST.hdf5"
    # infilename = f"/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrectionsByHelicity/PDFsFromCorrs/w_z_gen_dists_{corrname}_maxFiles_m1_pdfByHelicity_skimmed.hdf5"
    infilename = f"/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr_maxFiles_m1_msht20_msht20_pdfByHelicity.hdf5"
    
    

    outfolder = f"/work/submit/areimers/wmass/plots/genhistmaker/Z/ForAlphaS/WRemDev/NewSteer/variations/{corrname}"
    # outfolder = f"/work/submit/areimers/wmass/plots/genhistmaker/Z/ForAlphaS/WRemDev/NewSteer/variations/{corrname}_TEST"

    with h5py.File(infilename, "r") as h5file:
        results = input_tools.load_results_h5py(h5file)

    print(results.keys())

    samplename = "ZmumuPostVFP"

    print(results[samplename]["output"].keys())


    # histname = "nominal_gen_scetlib_dyturboMSHT20_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboCT18Z_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboMSHT20VarsCorr"
    histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr"
    hist = input_tools.read_and_scale(infilename, samplename, histname)
    print(hist)
    

    varaxisname = "vars"
    project_on = "ptVgen"

    # vars = [f"pdf0", "pdf2", "pdf5"]
    # vars = [0, 1, 2]
    # vars = [f"pdf{idx}" for idx in range(1, 59)]
    # vars = []
    # for idx in range(1, 30):
    #     vars += [f"pdf{idx}CT18ZDown", f"pdf{idx}CT18ZUp"]
    vars = [idx for idx in range(1, 65)]

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
                # rrange=[0.995, 1.005],
                nlegcols=1,
                xlim=[0, 100.], 
                ylim=[0, 90.], 
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