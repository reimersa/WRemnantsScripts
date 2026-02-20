from utilities.io_tools import input_tools

import h5py # type: ignore
import rabbit.io_tools # type: ignore




def main():



    filename = "/scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboCorr_OldNP_AllCT18Z_PDFFromCorr_ASFromCorr/fitresults_asimov.hdf5"

    fitresult, meta = rabbit.io_tools.get_fitresult(filename, "asimov", meta=True)

    print(fitresult["physics_models"]["Project ch0 ptll"]["channels"]["ch0"]["hist_prefit_inclusive_variations"].get())

    print(results["results_asimov"]["physics_models"].keys())


    # histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboMSHT20_pdfasCorr"
    # histname = "nominal_gen_pdfCT18Z"
    histname = "nominal_scetlib_dyturboCorr"
    hist = input_tools.read_and_scale(filename, samplename, histname)
    print(hist)
    print(hist.project("ptll"))
    

if __name__ == "__main__":
    main()