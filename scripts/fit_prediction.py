import argparse
import os
import pickle
import re

import lz4.frame

from utilities import common
from wremnants import theory_corrections
from wums import boostHistHelpers as hh
from wums import logging, output_tools, plot_tools
from scripts.corrections.make_theory_corr import read_corr
import matplotlib.pyplot as plt
from utilities.io_tools import input_tools
import mplhep
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
from iminuit import Minuit
from iminuit.util import describe
from functools import partial
from scipy.stats import chi2 as chi2_dist
from scipy.interpolate import make_smoothing_spline
import hist

def expxpow1_func(x, norm, pow, expo, const):
    return norm * np.power(x + const, pow) * np.exp(-expo*x)

def expxpow2_func(x, norm, pow, expo, expo2, const):
    return norm * np.power(x + const, pow) * np.exp(-expo*x + expo2*np.power(x, 2))

def expxpow3_func(x, norm, pow, expo, expo2, expo3, const):
    return norm * np.power(x + const, pow) * np.exp(-expo*x + expo2*np.power(x, 2) + expo3*np.power(x, 3))

unit_per_var = {
    "qT": " [GeV]",
    "Y": ""
}


def compute_uncertainty(func, x, minuit, eps=1e-7):
    """
    Given a model `func(x, *params)`, an array x of points,
    and a converged iminuit Minuit object, compute the 1sigma
    band via linear error propagation (finite-difference).
    
    Returns y_fit, y_unc arrays (both length len(x)).
    """
    # pull out parameter names, values & covariance
    names = minuit.parameters
    p0    = np.array([minuit.values[n] for n in names])
    cov   = np.array([[minuit.covariance[i][j] for j in names] for i in names])
    
    # central model
    y0 = func(x, *p0)
    
    # prepare gradient matrix (len(x) × npar)
    grads = np.zeros((len(x), len(p0)))
    
    # finite‐difference for each param
    for i in range(len(p0)):
        dp     = eps * max(abs(p0[i]), 1.0)      # step size
        p_up   = p0.copy();   p_up[i] += dp
        y_up   = func(x, *p_up)
        grads[:,i] = (y_up - y0) / dp
    
    # propagate: σ²(x) = grads(x) · cov · grads(x)^T
    var = np.einsum('ni,ij,nj->n', grads, cov, grads)
    return y0, np.sqrt(var)

def main():


    y_list_to_plot = [-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ]
    # y_list_to_plot = [1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ]
    # y_list_to_plot = [1.5,  1.75]
    # y_list_to_plot = []

    bins_ybins_postfixes_plottags = []
    # bins_ybins_postfixes_plottags += [({"vars" :  0}, (y_list_to_plot[0], y_list_to_plot[-1]), "Yinclusive", "(Y inclusive)")] 
    bins_ybins_postfixes_plottags += [
         ({"vars" :  0, "Y": ybin}, (y_list_to_plot[ybin], y_list_to_plot[ybin+1]), f"Y{y_list_to_plot[ybin]}_{y_list_to_plot[ybin+1]}".replace("-", "m").replace(".", "p"), f"(Y: [{y_list_to_plot[ybin]}, {y_list_to_plot[ybin+1]}])") for ybin in range(len(y_list_to_plot)-1)
    ]

    genpath = "/work/submit/areimers/wmass/TheoryCorrections"
    name_fo = "NNLOjet/Z/ZjNNLO/final/ptz"
    outfoldername = "/work/submit/areimers/wmass/plots/TheoryCorrections/Z/corrections_standalone/N3LOSmoothing"
    os.makedirs(outfoldername, exist_ok=True)
    qtmin = 7

    fitresults = []
    fitted_xvals = []
    y_edges = []
    for (bin, (yedgelow, yedgehigh), postfix, plottag) in bins_ybins_postfixes_plottags:
        print(f"--> Fitting {plottag}")
        var = "qT"
        y_edges.append((yedgelow, yedgehigh))

        y_list = [-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ]
        hist_fo = input_tools.read_nnlojet_pty_hist(os.path.join(genpath, name_fo), ybins=np.array(y_list), charge=0)
        # print(hist_fo[{"vars" :  0}])
        hist_fo_proj = hist_fo[bin].project(var)
        # hist_fo_proj = hist_fo[{"vars" :  0, "Y":0}].project(var)
        # hist_fo_proj = hist_fo[{"vars" : 0}].project(var)

        edges_orig   = hist_fo_proj.axes[0].edges
        counts_orig  = hist_fo_proj.values()
        errors_orig  = np.sqrt(hist_fo_proj.variances())

        rebin_factor = 1

        counts  = np.add.reduceat(counts_orig, np.arange(0, len(counts_orig), rebin_factor)) / rebin_factor
        errors    = np.sqrt(np.add.reduceat(errors_orig**2, np.arange(0, len(errors_orig), rebin_factor))) / rebin_factor
        edges   = edges_orig[::rebin_factor]
        centers = 0.5*(edges[:-1] + edges[1:])

        mask = (counts > 0) & (centers > qtmin)

        def chi2_expxpow1_minuit(norm, pow, expo, const):
            return np.sum(((counts[mask] - expxpow1_func(x=centers[mask], norm=norm, pow=pow, expo=expo, const=const)) / errors[mask])**2)

        def chi2_expxpow2_minuit(norm, pow, expo, expo2, const):
            return np.sum(((counts[mask] - expxpow2_func(x=centers[mask], norm=norm, pow=pow, expo=expo, expo2=expo2, const=const)) / errors[mask])**2)

        def chi2_expxpow3_minuit(norm, pow, expo, expo2, expo3, const):
            return np.sum(((counts[mask] - expxpow3_func(x=centers[mask], norm=norm, pow=pow, expo=expo, expo2=expo2, expo3=expo3, const=const)) / errors[mask])**2)

        initial_pars={
            "norm": 3e4,
            "pow": -2,
            "expo": 1e-2,
            "expo2": -3e-5,
            "expo3": -5e-8,
            "const": 10
        }
        m1 = Minuit(chi2_expxpow1_minuit, norm=initial_pars["norm"], pow=initial_pars["pow"], expo=initial_pars["expo"], const=initial_pars["const"])
        # m2 = Minuit(chi2_expxpow2_minuit, norm=initial_pars["norm"], pow=initial_pars["pow"], expo=initial_pars["expo"], expo2=initial_pars["expo2"], const=initial_pars["const"])
        # m3 = Minuit(chi2_expxpow3_minuit, norm=initial_pars["norm"], pow=initial_pars["pow"], expo=initial_pars["expo"], expo2=initial_pars["expo2"], expo3=initial_pars["expo3"], const=initial_pars["const"])

        fmin1 = m1.migrad(ncall=50000)
        m1.hesse()
        m1.minos()
        # m2.strategy = 1
        # m2.tol = 1e-5
        # fmin2 = m2.migrad(ncall=50000)
        # m2.hesse()
        # m2.minos()
        # fmin3 = m3.migrad(ncall=50000)
        # m3.hesse()
        # m3.minos()

        fitresults.append((expxpow1_func, m1))
        fitted_xvals.append(centers[mask])

        yfit1, y1_unc = compute_uncertainty(expxpow1_func, centers[mask], m1)
        # yfit2, y2_unc = compute_uncertainty(expxpow2_func, centers[mask], m2)
        # yfit3, y3_unc = compute_uncertainty(expxpow3_func, centers[mask], m3)


        print(fmin1)

        chi21 = m1.fval
        ndf1  = mask.sum() - m1.nfit
        # chi22 = m2.fval
        # ndf2  = mask.sum() - m2.nfit
        # chi23 = m3.fval
        # ndf3  = mask.sum() - m3.nfit

        fig = plt.figure(figsize=(6,6))
        gs  = fig.add_gridspec(2,1, height_ratios=[3,1], hspace=0.05)

        ax_main  = fig.add_subplot(gs[0])
        ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)


        datapoints = ax_main.errorbar(centers[mask], counts[mask], yerr=errors[mask], label="N3LO FO", fmt="o", ms=3, color="black")
        spline = make_smoothing_spline(centers[mask], counts[mask], lam=3.)
        spline_line, = ax_main.plot(centers[mask], spline(centers[mask]), '-', lw=2, label='spline')
        fit_line1, = ax_main.plot(centers[mask], yfit1, '-', lw=2, label=f"expxpow1 ({chi21:.1f}/{ndf1}, {chi2_dist.sf(chi21, ndf1)*100:.0f}%)")
        # fit_line2, = ax_main.plot(centers[mask], yfit2, '-', lw=2, label=f"expxpow2 ({chi22:.1f}/{ndf2}, {chi2_dist.sf(chi22, ndf2)*100:.0f}%)")
        # fit_line3, = ax_main.plot(centers[mask], yfit3, '-', lw=2, label=f"expxpow3 ({chi23:.1f}/{ndf3}, {chi2_dist.sf(chi23, ndf3)*100:.0f}%)")

        unc_scale = 1
        ax_main.fill_between(centers[mask], yfit1 - unc_scale*y1_unc, yfit1 + unc_scale*y1_unc, alpha=0.3, color=fit_line1.get_color())
        # ax_main.fill_between(centers[mask], yfit2 - unc_scale*y2_unc, yfit2 + unc_scale*y2_unc, alpha=0.3, color=fit_line2.get_color())
        # ax_main.fill_between(centers[mask], yfit3 - unc_scale*y3_unc, yfit3 + unc_scale*y3_unc, alpha=0.3, color=fit_line3.get_color())


        ax_main.set_ylabel("Prediction")
        ax_main.grid(True)
        handles, labels = ax_main.get_legend_handles_labels()
        ax_main.legend(handles[::-1], labels[::-1], loc="upper right", fontsize=14)
        # ax_main.legend(loc="upper right", fontsize=14)
        ax_main.set_yscale("log")

        # ratio = data / fit
        ratio_to = yfit1

        ratio     = counts[mask] / ratio_to
        ratio_err = errors[mask] / ratio_to
        ratio_fit1 = yfit1 / ratio_to
        ratio_fit1_unc = y1_unc / ratio_to
        # ratio_fit2 = yfit2 / ratio_to
        # ratio_fit2_unc = y2_unc / ratio_to
        # ratio_fit3 = yfit3 / ratio_to
        # ratio_fit3_unc = y3_unc / ratio_to
        ratio_spline = spline(centers[mask]) / ratio_to

        ax_ratio.errorbar(centers[mask], ratio, yerr=ratio_err, fmt="o", ms=3, color="black")
        ax_ratio.plot(centers[mask], ratio_spline, "-", lw=2, color=spline_line.get_color())
        fit_line_ratio1, = ax_ratio.plot(centers[mask], ratio_fit1, lw=2, color=fit_line1.get_color())
        # fit_line_ratio2, = ax_ratio.plot(centers[mask], ratio_fit2, lw=2)
        # fit_line_ratio3, = ax_ratio.plot(centers[mask], ratio_fit3, lw=2)
        ax_ratio.fill_between(centers[mask], ratio_fit1 - unc_scale*ratio_fit1_unc, ratio_fit1 + unc_scale*ratio_fit1_unc, alpha=0.3, color=fit_line_ratio1.get_color())
        # ax_ratio.fill_between(centers[mask], ratio_fit2 - unc_scale*ratio_fit2_unc, ratio_fit2 + unc_scale*ratio_fit2_unc, alpha=0.3, color=fit_line_ratio2.get_color())
        # ax_ratio.fill_between(centers[mask], ratio_fit3 - unc_scale*ratio_fit3_unc, ratio_fit3 + unc_scale*ratio_fit3_unc, alpha=0.3, color=fit_line_ratio3.get_color())
        ax_ratio.set_xlabel(f"{var}{unit_per_var[var]}")
        ax_ratio.set_ylabel("N3LO / fit")
        ax_ratio.set_ylim(0.90,1.10)
        ax_ratio.grid(True)

        # hide x-tick labels on main
        plt.setp(ax_main.get_xticklabels(), visible=False)


        # margins
        fig.subplots_adjust(
            left=0.17,
            right=0.98,
            top=0.98,
            bottom=0.15
        )
        plt.savefig(f"{outfoldername}/fit_N3LO_{var}_{postfix}.pdf")
        plt.close()


    # assemble new prediction for vs. y from N fits vs. qT
    smooth_hist = hist.Hist(
        hist.axis.Variable(y_list, name="Y"),  # 50 bins from –4 to +4
        hist.axis.Regular(100-qtmin, qtmin, 100, name="qT"),  # 50 bins from –4 to +6
        storage=hist.storage.Weight()            # optional: choose storage type
    )
    vars = np.zeros(smooth_hist.shape)

    for idx, (f, m) in enumerate(fitresults): # 1 fit per y bin
        x = fitted_xvals[idx]
        fitted_pars = np.array([m.values[n] for n in m.parameters])
        pred_smooth = f(x, *fitted_pars) # in this y bin
        y_center = (y_edges[idx][1] - y_edges[idx][0])/2 + y_edges[idx][0]
        y_arr = np.full_like(x, y_center)
        smooth_hist.fill(y_arr, x, weight=pred_smooth)


        _, errs = compute_uncertainty(f, x, m)
        vars[idx, :] = np.array(errs)**2 # 1d at this point, vs. x. 
        # smooth_hist.variances(flow=False)[...] = np.
    # print(smooth_hist.project("qT"))
    smooth_hist.variances(flow=False)[...] = vars

    y_list = [-5., -4., -3.5,  -3.25, -3., -2.75, -2.5,  -2.25 ,-2., -1.75, -1.5,  -1.25, -1., -0.75, -0.5,  -0.25,  0.,  0.25,  0.5, 0.75,  1.,  1.25,  1.5, 1.75, 2.,  2.25,  2.5, 2.75,  3.,  3.25,  3.5, 4.,  5.  ]
    hist_orig = input_tools.read_nnlojet_pty_hist(os.path.join(genpath, name_fo), ybins=np.array(y_list), charge=0)[{"vars": 0, "charge": 0, "qT": slice(qtmin, 100)}]
    sliced_hist_orig = hist.Hist(
        hist.axis.Variable(y_list, name="Y"),  # 50 bins from –4 to +4
        hist.axis.Regular(100-qtmin, qtmin, 100, name="qT"),  # 50 bins from –4 to +6
        storage=hist.storage.Weight()            # optional: choose storage type
    )
    sliced_hist_orig.values(flow=False)[...] = hist_orig.values(flow=False)
    sliced_hist_orig.variances(flow=False)[...] = hist_orig.variances(flow=False)
    print(sliced_hist_orig)


    colors = [f"{mcolors.to_hex(c)}" for c in list(plt.get_cmap("tab10").colors)]
    for var_to_compare in ["qT", "Y"]:
        fig = plot_tools.makePlotWithRatioToRef(
            hists=[smooth_hist.project(var_to_compare), sliced_hist_orig.project(var_to_compare)],
            labels=["Fitted", "Unsmoothed"],
            colors=colors[:2],
            xlabel=var_to_compare, 
            ylabel="Prediction",
            rlabel="unsm. / fit",
            rrange=[0.95, 1.05],
            nlegcols=1,
            xlim=None, 
            binwnorm=1.0, 
            baseline=True, 
            yerr=True,
            yerr_ratio=True,
            ratio_legend=False,
            linewidth=1,
        )
        fig.savefig(f"{outfoldername}/fitted_vs_unsmoothed_N3LO_{var_to_compare}.pdf")



    # 
    # hist[bh.underflow] = 0



if __name__ == "__main__":
    main()