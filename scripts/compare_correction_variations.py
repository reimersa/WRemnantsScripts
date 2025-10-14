import os
import re


from wums import plot_tools
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import numpy as np
import matplotlib.colors as mcolors


def corr_name(corrf):
    if not corrf.endswith(".pkl.lz4"):
        raise ValueError(f"File {corrf} is not a lz4 compressed pickle file")

    match = re.match(r"(.*)Corr[Z|W]\.pkl\.lz4", os.path.basename(corrf))
    return match[1]



# Define the exponential model: A * exp(-ax)
def exp_func(x, A, a):
    return A * np.exp(a*x)


def main():


    genpath = "/work/submit/areimers/wmass/TheoryCorrections"
    
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_higherprecision_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_variations/inclusive_Z_CT18Z_N3+0LL_frank_npvars_combined.pkl"
    # filename_resum = "SCETlib/ct18z_nplambda_n3+0ll_arnetest/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_higherprecision_combined.pkl"
    
    # filename_resum = "SCETlib/ct18z_nplambda_n3+0ll/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_reproducesold/inclusive_Z_CT18Z_N3+0LL_olduncertainties_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_reproducesold/inclusive_Z_CT18Z_N3+0LL_olduncertainties_Lambda24_higherprecision_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_reproducesold/inclusive_Z_CT18Z_N3+0LL_olduncertainties_cnu_omeganu_combined.pkl"

    # filename_resum = "SCETlib/ct18z_nplambda_n3+0ll/inclusive_Z_CT18Z_nplambda_N3+0LL_combined.pkl"
    # filename_resum = "SCETlib/ct18z_newnps_n3+0ll_reproducesold/inclusive_Z_CT18Z_N3+0LL_olduncertainties_final_allvars_higherprecision_combined.pkl"
    filename_resum = "SCETlib/ct18z_newnps_n3+0ll_lattice/inclusive_Z_CT18Z_N3+0LL_lattice_allvars_higherprecision_combined.pkl"

    variations_and_legnames = [
        # # Lambda2
        # ("pdf0", "nominal ($\\Lambda_{2} = 0.25$)"),
        # ("lambda20.0", "$\\Lambda_{2} = 0.0$"),
        # ("lambda20.5", "$\\Lambda_{2} = 0.5$"),
        # ("lambda2-0.25", "$\\Lambda_{2} = -0.25$"),
        # ("lambda20.75", "$\\Lambda_{2} = 0.75$"),
        # ("lambda2-0.5", "$\\Lambda_{2} = -0.50$"),
        # ("lambda21.0", "$\\Lambda_{2} = 1.00$"),

        # # Lambda4
        # ("pdf0", "nominal ($\\Lambda_{4} = 0.0625$)"),
        # ("lambda40.0", "$\\Lambda_{4} = 0.000$"),
        # ("lambda40.001", "$\\Lambda_{4} = 0.001$"),
        # ("lambda40.125", "$\\Lambda_{4} = 0.125$"),

        # # DeltaLambda2
        # ("pdf0", "nominal ($\\Delta\\Lambda_{2} = 0.125$)"),
        # ("delta_lambda20.105", "$\\Delta\\Lambda_{2} = 0.105$"),
        # ("delta_lambda20.145", "$\\Delta\\Lambda_{2} = 0.145$"),

        # lambda2_nu
        # ("pdf0", "nominal ($\\lambda_{2}^{\\nu} = 0.1$)"),
        # ("lambda2_nu0.0", "$\\lambda_{2}^{\\nu} = 0.0$"),
        # ("lambda2_nu0.2", "$\\lambda_{2}^{\\nu} = 0.2$"),
        # ("lambda2_nu-0.1", "$\\lambda_{2}^{\\nu} = -0.1$"),
        # ("lambda2_nu0.3", "$\\lambda_{2}^{\\nu} = 0.3$"),
        # ("lambda2_nu-0.2", "$\\lambda_{2}^{\\nu} = -0.2$"),
        # ("lambda2_nu0.4", "$\\lambda_{2}^{\\nu} = 0.4$"),
        # ("lambda2_nu-0.3", "$\\lambda_{2}^{\\nu} = -0.3$"),
        # ("lambda2_nu0.5", "$\\lambda_{2}^{\\nu} = 0.5$"),
        # ("lambda2_nu-0.4", "$\\lambda_{2}^{\\nu} = -0.4$"),
        # ("lambda2_nu0.6", "$\\lambda_{2}^{\\nu} = 0.6$"),

        # lambda4_nu
        # ("pdf0", "nominal ($\\lambda_{4}^{\\nu} = 0.01$)"),
        # ("lambda4_nu0.0", "$\\lambda_{4}^{\\nu} = 0.00$"),
        # # ("lambda4_nu0.0001", "$\\lambda_{4}^{\\nu} = 0.0001$"),
        # ("lambda4_nu0.02", "$\\lambda_{4}^{\\nu} = 0.02$"),
        # ("lambda4_nu0.03", "$\\lambda_{4}^{\\nu} = 0.03$"),


        ### LATTICE
        # individual variations
        # ("pdf0", "nominal (lattice)"),
        # ("lambda2_nu0.0538", "$\\lambda_{2}^{\\nu} = 0.0538$"),
        # ("lambda2_nu0.1202", "$\\lambda_{2}^{\\nu} = 0.1202$"),
        # ("lambda4_nu0.0008", "$\\lambda_{4}^{\\nu} = 0.0008$"),
        # ("lambda4_nu0.014", "$\\lambda_{4}^{\\nu} = 0.014$"),
        # ("lambda_inf_nu1.1784", "$\\lambda_{\\infty}^{\\nu} = 1.1784$"),
        # ("lambda_inf_nu2.1922", "$\\lambda_{\\infty}^{\\nu} = 2.1922$"),

        # eigenvariations
        # ("pdf0", "nominal (lattice)"),
        # ("lambda2_nu0.0696-lambda4_nu0.0122-lambda_inf_nu1.1Ext", "Eigenvar. 1 (up)"),
        # ("lambda2_nu0.1044-lambda4_nu0.0026-lambda_inf_nu2.1Ext", "Eigenvar. 1 (down)"),
        # ("lambda2_nu0.1153-lambda4_nu0.0032-lambda_inf_nu1.6Ext", "Eigenvar. 2 (up)"),
        # ("lambda2_nu0.0587-lambda4_nu0.0116-lambda_inf_nu1.6Ext", "Eigenvar. 2 (down)"),
        # ("lambda2_nu0.0873-lambda4_nu0.0092", "Eigenvar. 3 (up)"),
        # ("lambda2_nu0.0867-lambda4_nu0.0056", "Eigenvar. 3 (down)"),

        # lambda2 variations
        # ("pdf0", "nominal ($\\Lambda_{2} = 0.25$)"),
        # ("lambda20.0", "$\\Lambda_{2} = 0.0$"),
        # ("lambda20.5", "$\\Lambda_{2} = 0.5$"),
        # ("lambda2-0.25", "$\\Lambda_{2} = -0.25$"),
        # ("lambda20.75", "$\\Lambda_{2} = 0.75$"),
        # ("lambda2-0.5", "$\\Lambda_{2} = -0.5$"),
        # ("lambda21.0", "$\\Lambda_{2} = 1.0$"),

        # delta_lambda2 variations
        # ("pdf0", "nominal ($\\Delta\\Lambda_{2} = 0.125$)"),
        # ("delta_lambda20.085", "$\\Delta\\Lambda_{2} = 0.085$"),
        # ("delta_lambda20.165", "$\\Delta\\Lambda_{2} = 0.165$"),
        # ("delta_lambda20.065", "$\\Delta\\Lambda_{2} = 0.65$"),
        # ("delta_lambda20.185", "$\\Delta\\Lambda_{2} = 0.185$"),

        # lambda4 variations
        # ("pdf0", "nominal ($\\Lambda_{4} = 0.06$)"),
        # ("lambda40.04", "$\\Lambda_{4} = 0.04$"),
        # ("lambda40.08", "$\\Lambda_{4} = 0.08$"),
        # ("lambda40.0", "$\\Lambda_{4} = 0.0$"),
        # ("lambda40.12", "$\\Lambda_{4} = 0.12$"),
        # ("lambda4-0.06", "$\\Lambda_{4} = -0.06$"),
        # ("lambda40.18", "$\\Lambda_{4} = 0.18$"),
        # ("lambda4-0.18", "$\\Lambda_{4} = -0.18$"),
        # ("lambda40.3", "$\\Lambda_{4} = 0.3$"),

        # lambda_inf variations
        # ("pdf0", "nominal ($\\Lambda_{\\infty} = 1.0$)"),
        # ("lambda_inf0.9", "$\\Lambda_{\\infty} = 0.9$"),
        # ("lambda_inf1.1", "$\\Lambda_{\\infty} = 1.1$"),
        # ("lambda_inf0.8", "$\\Lambda_{\\infty} = 0.8$"),
        # ("lambda_inf1.2", "$\\Lambda_{\\infty} = 1.2$"),
        # ("lambda_inf0.5", "$\\Lambda_{\\infty} = 0.5$"),
        # ("lambda_inf1.5", "$\\Lambda_{\\infty} = 1.5$"),

        # new NP variations to reproduce old NP variations
        # ("pdf0", "nominal"),
        # ("lambda2-0.25", "$\\Lambda_{2}^{new} = -0.25$"),
        # ("lambda20.25", "$\\Lambda_{2}^{new} = 0.25$"),
        # ("lambda4.01", "$\\Lambda_{4}^{new} = 0.01$"),
        # ("lambda4.16", "$\\Lambda_{4}^{new} = 0.16$"),
        # ("delta_lambda2-0.02", "$\\Delta\\Lambda_{2}^{new} = -0.02$"),
        # ("delta_lambda20.02", "$\\Delta\\Lambda_{2}^{new} = 0.02$"),



        # old NP variations (target for new NP model)
        # ("pdf0", "nominal"),
        
        # # cnu + omeganu part
        # ("omega_nu0.5", "$(c_{\\nu}, \\omega_{\\nu}) = (nom. (0.1), 0.5)$"),
        # ("c_nu-0.1-omega_nu0.5", "$(c_{\\nu}, \\omega_{\\nu}) = (-0.1, 0.5)$"),

        # # lambda2/4 + dLambda2
        # ("Lambda20.25", "$\\Lambda_{2}^{old} = 0.25$"),
        # ("Lambda2-0.25", "$\\Lambda_{2}^{old} = -0.25$"),
        # ("Lambda4.01", "$\\Lambda_{4}^{old} = 0.01$"),
        # ("Lambda4.16", "$\\Lambda_{4}^{old} = 0.16$"),
        # ("Delta_Lambda20.02", "$\\Delta\\Lambda_{2}^{old} = 0.02$"),
        # ("Delta_Lambda2-0.02", "$\\Delta\\Lambda_{2}^{old} = -0.02$"),



        # # new NP variations to reproduce old NP variations
        # ("pdf0", "nominal"),

        # # # cnu + omeganu part
        # ("lambda2_nu-0.5", "$(\\lambda_{2, \\nu}^{new}, \\lambda_{\\infty, \\nu}^{new}) = (-0.5, nom. (-0.2))$"),
        # ("lambda2_nu0.5-lambda_inf_nu0.2", "$(\\lambda_{2, \\nu}^{new}, \\lambda_{\\infty, \\nu}^{new}) = (0.5, 0.2)$"),

        # # # lambda2/4 + dLambda2
        # ("lambda2-0.25", "$\\Lambda_{2}^{new} = -0.25$"),
        # ("lambda20.25", "$\\Lambda_{2}^{new} = 0.25$"),
        # ("lambda4.01", "$\\Lambda_{4}^{new} = 0.01$"),
        # ("lambda4.16", "$\\Lambda_{4}^{new} = 0.16$"),
        # ("delta_lambda2-0.02", "$\\Delta\\Lambda_{2}^{new} = -0.02$"),
        # ("delta_lambda20.02", "$\\Delta\\Lambda_{2}^{new} = 0.02$"),


        # new NP with lattice central values + variations
        ("pdf0", "nominal (lattice)"),
        # ("lambda2_nu0.0696-lambda4_nu0.0122-lambda_inf_nu1.1Ext", "Eigenvar. 1 (up)"),
        # ("lambda2_nu0.1044-lambda4_nu0.0026-lambda_inf_nu2.1Ext", "Eigenvar. 1 (down)"),
        # ("lambda2_nu0.1153-lambda4_nu0.0032-lambda_inf_nu1.6Ext", "Eigenvar. 2 (up)"),
        # ("lambda2_nu0.0587-lambda4_nu0.0116-lambda_inf_nu1.6Ext", "Eigenvar. 2 (down)"),
        # ("lambda2_nu0.0873-lambda4_nu0.0092", "Eigenvar. 3 (up)"),
        # ("lambda2_nu0.0867-lambda4_nu0.0056", "Eigenvar. 3 (down)"),

        ("lambda20.0", "$\\Lambda_{2}^{new} = 0.00$"),
        ("lambda20.5", "$\\Lambda_{2}^{new} = 0.50$"),
        ("lambda40.01", "$\\Lambda_{4}^{new} = 0.01$"),
        ("lambda40.16", "$\\Lambda_{4}^{new} = 0.16$"),
        ("delta_lambda20.105", "$\\Delta\\Lambda_{2}^{new} = 0.105$"),
        ("delta_lambda20.145", "$\\Delta\\Lambda_{2}^{new} = 0.145$"),



    ]

    print(input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum)))

    hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum))[{"vars" : var}] for (var, _) in variations_and_legnames]
    # hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum))[{"vars" : var, "qT": 0}] for (var, _) in variations_and_legnames]
    # hists_resum = [input_tools.read_scetlib_hist(os.path.join(genpath, filename_resum))[{"vars" : var, "Y": 3}] for (var, _) in variations_and_legnames]
    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone"
    os.makedirs(outfoldername, exist_ok=True)
    # print(hists_resum[0][{"Y": 0}].project("qT"))


    variables = ["qT", "Y"]
    # variables = ["qT"]
    # variables = ["Y"]

    for v in variables:
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[hist_resum.project(v) for hist_resum in hists_resum],
            labels=[vl[1] for vl in variations_and_legnames],
            colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab20").colors)][1:len(variations_and_legnames)+1],
            xlabel=f"{v} (GeV)" if v=="qT" else v, 
            ylabel="Resumm. prediction",
            rlabel=f"x/nominal",
            rrange=[0.95, 1.05] if v=="qT" else [0.9990, 1.0010],
            nlegcols=1,
            xlim=[0., 100.] if v == "qT" else [-5., 5.], 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
            logy=False,
        )
        fig.savefig(f"{outfoldername}/correction_comparison_variations_resum_{v}.pdf")


        

if __name__ == "__main__":
    main()