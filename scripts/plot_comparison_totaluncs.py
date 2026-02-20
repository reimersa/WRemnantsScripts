import rabbit.io_tools # type: ignore
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.ticker import AutoMinorLocator # type: ignore




def main():

    class setting:
        def __init__(self, filename, legname, color):
            self.filename = filename
            self.legname = legname
            self.color = color


    varname = "ptll"
    xaxis_title = "$p_{\\mathrm{T}}^{\\mu\\mu}$ (GeV)"
    xlim=(0, 100)

    # varname = "yll"
    # xaxis_title = "$y^{\\mu\\mu}$"
    # xlim = (-5., 5.)

    settings = [
        setting(filename="/scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboCorr_OldNP_AllCT18Z_PDFFromCorr_ASFromCorr/fitresults_asimov.hdf5", legname="Old NP model (N$^{3{+}0}$LL+NNLO)", color="tab:gray"),
        setting(filename="/scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints/fitresults_asimov.hdf5", legname="New NP model: nominal (N$^{3{+}0}$LL+NNLO)", color="tab:red"), 
        setting(filename="/scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_NoLatticeConstraintsx10/fitresults_asimov.hdf5", legname="New NP model: alternative (N$^{3{+}0}$LL+NNLO)", color="tab:blue"),
    ]

    fitresults = [rabbit.io_tools.get_fitresult(x.filename, "asimov", meta=False) for x in settings]
    hists = [x["physics_models"][f"Project ch0 {varname}"]["channels"]["ch0"]["hist_prefit_inclusive"].get() for x in fitresults]
    axes = [h.axes for h in hists]
    reluncs = [h.variances() ** 0.5 / h.values() for h in hists]
    hist_as_vars = [fitresults[-1]["physics_models"][f"Project ch0 {varname}"]["channels"]["ch0"]["hist_prefit_inclusive_variations"].get()[{"vars": "pdfAlphaS", "downUpVar": x}] for x in [0, 1]]
    reluncs_as_vars = [x.values() / hists[-1].values() for x in hist_as_vars]


    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for idx, h in enumerate(hists):
        edges = np.asarray(h.axes[0] .edges) # bin edges length N+1
        y_lo = 0.0 - reluncs[idx]
        y_hi = 0.0 + reluncs[idx]
        x_step, lo_step = step_expand(edges, y_lo)
        _,      hi_step = step_expand(edges, y_hi)
        ax.fill_between(x_step, lo_step, hi_step, alpha=0.25, label=settings[idx].legname, color=settings[idx].color)
        ax.step(x_step, lo_step, where="mid", linewidth=1.5, color=settings[idx].color)
        ax.step(x_step, hi_step, where="mid", linewidth=1.5, color=settings[idx].color)
    edges_asvars = np.asarray(hist_as_vars[0].axes[0] .edges) # bin edges length N+1
    x_step, lo_step = step_expand(edges_asvars, reluncs_as_vars[0] - 1.0)
    _,      hi_step = step_expand(edges_asvars, reluncs_as_vars[1] - 1.0)
    print(x_step)
    print(lo_step)
    print(hi_step)
    # ax.fill_between(x_step, lo_step, hi_step, alpha=0.25, label="alphas", color="black")
    ax.step(x_step, lo_step, where="mid", linewidth=1.5, color="black", label="$\\alpha_{s}\\pm 0.002$")
    ax.step(x_step, hi_step, where="mid", linewidth=1.5, color="black")

    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(axis="both", which="major", labelsize=14, length=5, width=1.0)
    ax.tick_params(axis="both", which="minor", length=3, width=0.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))  # 4 minors between majors
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.axhline(0.0, linewidth=1, color="black", linestyle="--")
    ax.set_xlabel(xaxis_title, ha="right", fontsize=16)
    ax.xaxis.set_label_coords(1.0, -0.08)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylabel("Relative total uncertainty", ha="right", fontsize=16)
    ax.yaxis.set_label_coords(-0.095, 1.0)
    ax.legend(frameon=False, fontsize=12)

    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.14, top=0.96)

    outname = f"/work/submit/areimers/wmass/plots/fitresults/ForAlphaS/WRemDev/NewSteer/comparisons/comparison_totaluncs_{varname}.pdf"
    fig.savefig(outname)
    print(f"--> Wrote plot {outname}")




def step_expand(edges, y):
    """
    Convert bin edges (N+1) and per-bin y (N) into arrays suitable for stepwise plotting.
    """
    edges = np.asarray(edges)
    y = np.asarray(y)
    assert edges.ndim == 1 and y.ndim == 1
    assert len(edges) == len(y) + 1

    x_step = np.repeat(edges, 2)[1:-1]  # length 2N
    y_step = np.repeat(y, 2)            # length 2N
    return x_step, y_step
    

if __name__ == "__main__":
    main()