import os
import h5py

from wums import boostHistHelpers as hh
from wums import output_tools, plot_tools
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import mplhep
import numpy as np
import tensorflow as tf
from scipy.optimize import curve_fit
from wums.ioutils import pickle_load_h5py


def print_tree(name, obj):
    indent = "  " * name.count('/')
    if isinstance(obj, h5py.Dataset):
        print(f"{indent}{name}  ⟶  shape={obj.shape}, dtype={obj.dtype}")
    else:
        print(f"{indent}{name}/")

def main():
    plot_variations = True

    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_nominal/ZMassDilepton.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll__perBinStatCorrUnc/ZMassDilepton.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_perBinStatCorrUncNdim/ZMassDilepton.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_n4p0lln3lo_unsmoothed_n3pllfixed_TESTCORR/ZMassDilepton.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassWLike_eta_pt_charge_scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorr/ZMassWLike.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassWLike_eta_pt_charge_scetlib_nnlojet_N4p0LLN3LOUnsmoothed_N3pLLFixedCorr_perbinstatcorrunc/ZMassWLike.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassWLike_eta_pt_charge_scetlib_nnlojetN4p0LLN3LOCorr_SmoothFO_Kenneth_perbinstatcorrunc/ZMassWLike.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/WMass_charge_pt_eta_scetlib_nnlojet_N3p1LLN3LOUnsmoothed_N3pLLFixedCorr_FixedAlphaS_PerBinStatCorrUnc/WMass.hdf5"
    # infilename = f"/scratch/submit/cms/areimers/wmass/fitinputs/ZMassWLike_eta_pt_charge_scetlib_dyturboCorr_FixedAlphaS/ZMassWLike.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ZMassWLike_eta_pt_charge_scetlib_dyturboCorr_FixedAlphaS_BinnedScale//ZMassWLike.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_scetlib_dyturbo_NewNPModel_OldValsAndVarsCorr_FixedAlphaS_PerBinStatCorrUnc//ZMassDilepton.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_scetlib_dyturbo_NewNPModel_LatticeValsAndVarsCorr_FixedAlphaS_PerBinStatCorrUnc//ZMassDilepton.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ZMassDilepton_ptll_yll_scetlib_dyturboCorr_FixedAlphaS_AlphaS_PerBinStatCorrUnc//ZMassDilepton.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ForAlphaS/WRemDev/ZMassDilepton_ptll_yll_scetlib_dyturbo_NewNPModel_LatticeValsOldVarsCorr_FixedAlphaS_AlphaS//ZMassDilepton.hdf5"
    # infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ForAlphaS/WRemDev/ZMassDilepton_ptll_yll_scetlib_dyturboN3p0LL_LatticeNPCorr_FixedAlphaS_HackedPDFVarsFromHels_OwnPDFVars_AlphaS//ZMassDilepton.hdf5"
    infilename = "/scratch/submit/cms/areimers/wmass/fitinputs/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints//ZMassDilepton.hdf5"


    # outfolder = f"/work/submit/areimers/wmass/plots/fitinputs/ZMassDilepton_ptll_nominal"
    # outfolder = f"/work/submit/areimers/wmass/plots/fitinputs/ZMassDilepton_ptll__perBinStatCorrUnc"
    # outfolder = f"/work/submit/areimers/wmass/plots/fitinputs/ZMassDilepton_ptll_REPEAT"




    f = h5py.File(infilename, mode="r")
    f.visititems(print_tree)

    constraintweights = f['hconstraintweights']
    data_obs = f['hdata_obs'][...]
    logk = f['hlogk'][...]
    noigroupidxs = f['hnoigroupidxs'][...]
    noigroups = f['hnoigroups'][...]
    norm = f['hnorm'][...]
    procs = f['hprocs'][...]
    pseudodata = f['hpseudodata'][...]
    pseudodatanames = f['hpseudodatanames'][...]
    signals = f['hsignals'][...]
    sumw = f['hsumw'][...]
    sumw2 = f['hsumw2'][...]
    systgroupidxs = f['hsystgroupidxs'][...]
    systgroups = f['hsystgroups'][...]
    systs = f['hsysts'][...]
    systsnoconstraint = f['hsystsnoconstraint'][...]
    # systsnoprofile = f['hsystsnoprofile'][...]
    # meta = f['meta']
    meta = pickle_load_h5py(f["meta"])
    chan = meta["channel_info"]["ch0"]
    var = chan["axes"][0]

    # print("top-level keys:", list(f.keys()))
    # print("constraintweights", len(constraintweights), type(constraintweights), constraintweights.shape, constraintweights)
    # print("data_obs", len(data_obs), type(data_obs), data_obs.shape, data_obs)
    # print("logk", len(logk), type(logk), logk.shape, logk)
    # print("noigroupidxs", len(noigroupidxs), type(noigroupidxs), noigroupidxs.shape, noigroupidxs)
    # print("noigroups", len(noigroups), type(noigroups), noigroups.shape, noigroups)
    # print("norm", len(norm), type(norm), norm.shape, norm)
    # print("procs", len(procs), type(procs), procs.shape, procs)
    # print("pseudodata", len(pseudodata), type(pseudodata), pseudodata.shape, pseudodata)
    # print("pseudodatanames", len(pseudodatanames), type(pseudodatanames), pseudodatanames.shape, pseudodatanames)
    # print("signals", len(signals), type(signals), signals.shape, signals)
    # print("sumw", len(sumw), type(sumw), sumw.shape, sumw)
    # print("sumw2", len(sumw2), type(sumw2), sumw2.shape, sumw2)
    # print("systgroupidxs", len(systgroupidxs), type(systgroupidxs), systgroupidxs.shape, systgroupidxs)
    # print("systgroups", len(systgroups), type(systgroups), systgroups.shape, systgroups)
    print("systs", len(systs), type(systs), systs.shape, systs)
    # print("systsnoconstraint", len(systsnoconstraint), type(systsnoconstraint), systsnoconstraint.shape, systsnoconstraint)
    # print("systsnoprofile", len(systsnoprofile), type(systsnoprofile), systsnoprofile.shape, systsnoprofile)
    # # print("meta", len(meta), type(meta), meta)
    # # print(chan)
    
    print(var)
    print(len(systs))
    for s in systs:
        # if "theoryCorrStat" in str(s) or "per_bin_stat" in str(s):
        # if "scetlib" in str(s):
        # if "resum" in str(s):
        if "scale" in str(s).lower():
        # if "pdf" in str(s):
            print(f"--> Found relevant {str(s)} in list of systs")










    # compare statistical uncertainties
    # print(hist[{"vars": 0}])
    # relstatunc_nnlojet = hist[{"vars": 0}].sum().variance / hist[{"vars": 0}].sum().value
    # print(relstatunc_nnlojet)
    # sqrt(h.sum().variance)/h.sum().value

    
    

if __name__ == "__main__":
    main()