import os
import re


from wums import plot_tools
import matplotlib.pyplot as plt # type: ignore
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors # type: ignore

import lz4.frame # type: ignore
import pickle


def main():


    # correction = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorrZ.pkl.lz4"
    correction = "/work/submit/areimers/wmass/WRemnants/wremnants-data/data/TheoryCorrections/scetlib_dyturboMSHT20_pdfasCorrZ.pkl.lz4"

    variations_num_den_and_legnames = [
        # # pdfas vars for MSTH20 (called like normal PDF vars)
        ("pdf0", "as0118", "nominal"),
        ("pdf2", "as0116", "as down"),
        ("pdf5", "as0120", "as up"),
    ]

    # histname_numh = "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas_hist"
    histname_numh = "scetlib_dyturboMSHT20_pdfas_hist"







    hists_num = [pickle.load(lz4.frame.open(correction))["Z"][histname_numh][{"vars": var_num}] for var_num, _, _ in variations_num_den_and_legnames]
    histname_den = "minnlo_ref_hist"
    hists_den = [pickle.load(lz4.frame.open(correction))["Z"][histname_den][{"vars": var_den}] for _, var_den, _ in variations_num_den_and_legnames]
    print(hists_num)
    print(hists_den)

    variables = ["qT"]

    for v in variables:
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[h.project(v) for h in hists_num],
            labels=[l[2] for l in variations_num_den_and_legnames],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab20").colors)][1:len(variations_num_den_and_legnames)+1],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="Numerator",
            rlabel=f"x/nominal",
            rrange=[0.95, 1.05] if v=="qT" else [0.9990, 1.0010],
            nlegcols=1,
            xlim=[0., 100.] if v == "qT" else [-5., 5.], 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            # yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=False,
        )
        fig.savefig(f"/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone/correction_full_numerator_{v}.pdf")

        fig = plot_tools.makePlotWithRatioToRef(
            hists=[h.project(v) for h in hists_den],
            labels=[l[2] for l in variations_num_den_and_legnames],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab20").colors)][1:len(variations_num_den_and_legnames)+1],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="Denominator",
            rlabel=f"x/nominal",
            rrange=[0.95, 1.05] if v=="qT" else [0.9990, 1.0010],
            nlegcols=1,
            xlim=[0., 100.] if v == "qT" else [-5., 5.], 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            # yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=False,
        )
        fig.savefig(f"/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone/correction_full_denominator_{v}.pdf")


        

if __name__ == "__main__":
    main()