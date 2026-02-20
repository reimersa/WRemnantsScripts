import pickle
import lz4.frame




def main():


    # correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorrZ.pkl.lz4"
    # correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorrZ.pkl.lz4"
    # correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboN3p0LL_LatticeNP_MSHT20CorrZ.pkl.lz4"
    # correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboMSHT20_pdfasCorrZ.pkl.lz4"
    # correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboCT18Z_pdfasCorrZ.pkl.lz4"
    correctionname = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_CorrZ.pkl.lz4"
    
    corrfile = pickle.load(lz4.frame.open(correctionname))
    proc = "Z"
    print(corrfile[proc])
    print(corrfile[proc]["scetlib_dyturboCT18Z_pdfas_hist"])


if __name__ == "__main__":
    main()