from utilities.io_tools import input_tools

import h5py # type: ignore




def main():


    # filename = "/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_maxFiles_m1.hdf5"
    # filename = "/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_maxFiles_m1_ct18z_msht20.hdf5"
    # filename = "/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_maxFiles_m1_msht20_msht20.hdf5"
    # filename = "/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/NewSteer/mz_dilepton_scetlib_dyturboCorr_maxFiles_m1_OldNP_AllCT18Z.hdf5"
    # filename = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrectionsByHelicity/AlphaS/w_z_gen_dists_scetlib_dyturboMSHT20_pdfasCorr_maxFiles_m1_skimmed.hdf5"
    # filename = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrectionsByHelicity/AlphaS/w_z_gen_dists_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr_maxFiles_m1_skimmed.hdf5"
    # filename = "/scratch/submit/cms/areimers/wmass/histmaker/ForAlphaS/WRemDev/NewSteer/mz_dilepton_scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_Corr_maxFiles_m1_NewNP_FineBins_AllCT18Z.hdf5"
    filename = "/scratch/submit/cms/areimers/wmass/gendistributions/w_z_gen_dists_maxFiles_m1__finePtAbsyQ.hdf5"
    

    with h5py.File(filename, "r") as h5file:
        results = input_tools.load_results_h5py(h5file)

    print(results.keys())

    samplename = "ZmumuPostVFP"
    print(results[samplename]["output"].keys())


    # histname = "nominal_gen_scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr"
    # histname = "nominal_gen_scetlib_dyturboMSHT20_pdfasCorr"
    # histname = "nominal_gen_pdfCT18Z"
    # histname = "nominal_scetlib_dyturboCorr"
    histname = "gen_massVgen_scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_Corr"

    hist = input_tools.read_and_scale(filename, samplename, histname)
    print(hist)
    print(hist.project("ptll"))
    

if __name__ == "__main__":
    main()