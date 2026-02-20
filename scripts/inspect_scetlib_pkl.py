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
    plot_variations = False

    # filename = "/work/submit/areimers/wmass/TheoryCorrections/SCETlib/msht20nnlo_nplambda_pdfvars/inclusive_Z_MSHT20nnlo_nplambda_pdfas_nnlo_sing_combined.pkl"
    filename = "/work/submit/areimers/wmass/TheoryCorrections/SCETlib/ct18z_nptan2_n3p0ll/inclusive_Z_CT18Z_N3p0LL_newnp_frankvals.pkl"
    hist = input_tools.read_scetlib_hist(filename)

    print(hist)
    
    

if __name__ == "__main__":
    main()