import os
from dataclasses import dataclass, field
from enum import Enum









def main():


    outfolder_base = "/scratch/submit/cms/areimers/wmass"
    plotfolder_base = "/work/submit/areimers/wmass/plots"
    wremnants_folder = "/work/submit/areimers/wmass/WRemnants"

    max_threads = 350

    configs = [

        ##### New NP with lattice-constrained eigenvariations
        # {
        #     "theory_corrs": ["scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO", "scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_pdfas", "scetlib_dyturbo_LatticeNP_CT18Z_N3p0LL_N2LO_pdfvars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "WithLatticeConstraints",
        #     "np_model": "LatticeEigvars"
        # },
        {
            "theory_corrs": ["scetlib_dyturbo_LatticeNP_CT18Z_N2p1LL_N2LO", "scetlib_dyturbo_LatticeNP_CT18Z_N2p1LL_N2LO_pdfas", "scetlib_dyturbo_LatticeNP_CT18Z_N2p1LL_N2LO_pdfvars"],
            "pdfs": ["ct18z"],
            "postfix_histmaker": "",
            "pdf_inflation_factor": 1.0,
            "postfix_rabbit": "WithLatticeConstraints",
            "np_model": "LatticeEigvars"
        },
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP_MSHT20", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20Vars", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "NewNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.7,
        #     "postfix_rabbit": "NewNP_AllMSHT20_PDFFromCorr_ASFromCorr_PDFInflated1p7_WithLatticeConstraints",
        #     "np_model": "LatticeEigvars"
        # },

        # ##### New NP, no lattice constraints, instead inflating the individual vars x10
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP", "scetlib_dyturboN3p0LL_LatticeNP_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVars", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_NoLatticeConstraintsx10",
        #     "np_model": "LatticeNoConstraintsx10"
        # },
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP_MSHT20", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20Vars", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "NewNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.70,
        #     "postfix_rabbit": "NewNP_AllMSHT20_PDFFromCorr_ASFromCorr_PDFInflated1p7_NoLatticeConstraintsx10",
        #     "np_model": "LatticeNoConstraintsx10"
        # },
    ]



    # modes = [Mode.as_dilepton_1d]
    # modes = [Mode.as_dilepton_2d]
    modes = [Mode.as_dilepton_4d]




    for mode in modes:
        for config in configs:
            runner = AnalysisRunner(mode=mode, max_threads=max_threads, outfolder_base=outfolder_base, plotfolder_base=plotfolder_base, wremnants_folder=wremnants_folder, theory_corrs=config["theory_corrs"], pdfs=config["pdfs"], pdf_inflation_factor=config["pdf_inflation_factor"], np_model=config["np_model"], postfix_histmaker=config["postfix_histmaker"], postfix_rabbit=config["postfix_rabbit"])

            runner.run_histmaker()

            # runner.run_setup_rabbit()

            # runner.run_fit_rabbit()

            # runner.plot_fit_rabbit(prefit=True)
            # runner.plot_fit_rabbit(prefit=False) # postfit

            # runner.plot_impacts()

            # runner.plot_param_corr(theory_corr_from_hels=True)
            # runner.plot_param_corr(theory_corr_from_hels=True, np_model="LatticeNoConstraints")





class Mode(Enum):
    as_dilepton_1d = "as_dilepton_1d"
    as_dilepton_2d = "as_dilepton_2d"
    as_dilepton_4d = "as_dilepton_4d"

mode_configs = {
    Mode.as_dilepton_2d: {
        "exec_histmaker": "mz_dilepton.py",
        "args_histmaker": ["--axes yll ptll", "--csVarsHist"],
        "outputfile_histmaker": "mz_dilepton",
        "subfolder_histmaker": "AlphaS/AllCorrs",
        "args_setup_rabbit": ["--fitvar ptll-yll", "--noi alphaS"],
        "fitinput_folder": "ZMassDilepton_ptll_yll",
        "fitinput_filename": "ZMassDilepton.hdf5",
        "args_unblind": "",
    },
}
mode_configs[Mode.as_dilepton_4d] = mode_configs[Mode.as_dilepton_2d].copy()
mode_configs[Mode.as_dilepton_4d]["args_setup_rabbit"] = ["--fitvar ptll-yll-cosThetaStarll_quantile-phiStarll_quantile", "--noi alphaS"]
mode_configs[Mode.as_dilepton_4d]["fitinput_folder"] = "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile"

mode_configs[Mode.as_dilepton_1d] = mode_configs[Mode.as_dilepton_2d].copy()
mode_configs[Mode.as_dilepton_1d]["args_setup_rabbit"] = ["--fitvar ptll", "--noi alphaS"]
mode_configs[Mode.as_dilepton_1d]["fitinput_folder"] = "ZMassDilepton_ptll"


def make_postfix_plot(nominal_corr, postfix):
    result = f"{nominal_corr}Corr"
    if postfix:
        result += f"_{postfix}"
    return result


@dataclass
class AnalysisRunner:
    mode: Mode
    outfolder_base: str
    plotfolder_base: str
    wremnants_folder: str
    theory_corrs: list[str]
    pdf_inflation_factor: float
    np_model: str
    postfix_histmaker: str
    postfix_rabbit: str
    pdfs: list[str] = field(init=["ct18z"])
    max_threads: int = -1
    
    nominal_corr: str = field(init=False)
    postfix_fit_plot: str = field(init=False)
    exec_histmaker: str = field(init=False)
    args_histmaker: list[str] = field(init=False)
    outputfile_histmaker: str = field(init=False)
    args_setup_rabbit: list[str] = field(init=False)
    fitinput_folder: str = field(init=False)
    fitinput_filename: str = field(init=False)
    args_unblind: str = field(init=False)

    def __post_init__(self):
        self.nominal_corr = self.theory_corrs[0]

        self.postfix_histmaker = f"{self.postfix_histmaker}"
        if self.pdfs[0] != "ct18z":
            self.postfix_histmaker = f"{self.pdfs[0]}_{self.postfix_histmaker}"
        self.postfix_fit_plot = make_postfix_plot(nominal_corr=self.nominal_corr, postfix=self.postfix_rabbit)

        self.pdfs += ["msht20mcrange_renorm", "msht20mbrange_renorm"]

        self.exec_histmaker = mode_configs[self.mode]["exec_histmaker"]
        self.args_histmaker = mode_configs[self.mode]["args_histmaker"]
        self.outputfile_histmaker = mode_configs[self.mode]["outputfile_histmaker"]
        self.subfolder_histmaker = mode_configs[self.mode]["subfolder_histmaker"]
        self.args_setup_rabbit = mode_configs[self.mode]["args_setup_rabbit"]
        self.fitinput_folder = mode_configs[self.mode]["fitinput_folder"]
        self.fitinput_filename = mode_configs[self.mode]["fitinput_filename"]
        self.args_unblind = mode_configs[self.mode]["args_unblind"]

        self.outfolder_histmaker = os.path.join(self.outfolder_base, 'histmaker')
        if self.subfolder_histmaker:
            self.outfolder_histmaker = os.path.join(self.outfolder_histmaker, self.subfolder_histmaker)

        self.outfolder_setuprabbit = os.path.join(self.outfolder_base, 'fitinput')
        if self.subfolder_histmaker:
            self.outfolder_setuprabbit = os.path.join(self.outfolder_setuprabbit, self.subfolder_histmaker)

        self.outfolder_fit = os.path.join(self.outfolder_base, 'fitresults')
        if self.subfolder_histmaker:
            self.outfolder_fit = os.path.join(self.outfolder_fit, self.subfolder_histmaker)

        self.outfolder_plots = os.path.join(self.plotfolder_base, 'fitresults')
        if self.subfolder_histmaker:
            self.outfolder_plots = os.path.join(self.outfolder_plots, self.subfolder_histmaker)
        







        


    def run_histmaker(self):

        args = self.args_histmaker
        args.append(f"-o '{self.outfolder_histmaker}'")
        args.append("--maxFiles '-1'")
        args.append(f"-j {self.max_threads}" if self.max_threads > 0 else "")
        args.append("--theoryCorr " + " ".join(self.theory_corrs))
        args.append("--pdfs " + " ".join(self.pdfs))
        args.append(f"-p {self.postfix_histmaker}" if self.postfix_histmaker else "")

        args = [a for a in args if a]
        all_args = " ".join(args)
            
        command = f"python scripts/histmakers/{self.exec_histmaker} {all_args}"

        print(command)
        os.system(command)



    def run_setup_rabbit(self):
        args = self.args_setup_rabbit
        args.append("--verbose 3")
        args.append(f"-i '{os.path.join(self.outfolder_histmaker, f'{self.outputfile_histmaker}_{self.nominal_corr}_Corr_maxFiles_m1' + (f'_{self.postfix_histmaker}' if self.postfix_histmaker else '') + '.hdf5')}'")
        args.append(f"--npUnc '{self.np_model}'")
        args.append("--pdfUncFromCorr")
        args.append(f"-o '{self.outfolder_setuprabbit}'")
        args.append("--realData")
        args.append(f"--scalePdf {self.pdf_inflation_factor}")
        args.append(f"-p {self.postfix_fit_plot}")

        args = [a for a in args if a]
        all_args = " ".join(args)

        command = f"python scripts/rabbit/setupRabbit.py {all_args}"

        print(command)
        os.system(command)

    def run_fit_rabbit(self):
        args = [f"'{os.path.join(self.outfolder_setuprabbit, f'{self.fitinput_folder}_{self.postfix_fit_plot}', self.fitinput_filename)}'"]
        args.append("--globalImpacts")
        args.append("--computeHistErrors")
        args.append("--computeVariations")
        args.append("--saveHists")
        args.append("--saveHistsPerProcess")
        args.append("--computeHistCov")
        args.append("-m Project ch0 ptll")
        args.append("-m Project ch0 yll")
        args.append("-m Project ch0 cosThetaStarll_quantile")
        args.append("-m Project ch0 phiStarll_quantile")
        args.append("--doImpacts")
        args.append("-t -1 --postfix asimov")
        # args.append("-t 0 --postfix data")
        args.append(f"-o '{os.path.join(self.outfolder_fit, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'")
        args.append(self.args_unblind)

        args = [a for a in args if a]
        all_args = " ".join(args)


        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py {all_args}"
        print(command)
        os.system(command)


    def plot_fit_rabbit(self, prefit=True, show_data=False):
        args = [f"'{os.path.join(self.outfolder_fit, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}'"]
        args.append(f"-o '{os.path.join(self.outfolder_plots, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'")
        args.append("--config utilities/styles/styles.py")
        args.append("--title CMS")
        args.append("--titlePos 1")
        args.append("--subtitle WiP")
        args.append("--extraTextLoc 0.5 0.94")
        args.append("--showVariations both")
        args.append("--legCols 1")
        args.append("--yscale 1.2")
        args.append("--legCol 1")
        args.append("--legSize small")

        # args.append("--rrange 0.98 1.02")
        args.append("--rrange 0.90 1.10")

        args.append("-m Project ch0 ptll")
        # args.append("-m Project ch0 yll")
        # args.append("-m Project ch0 cosThetaStarll_quantile")
        # args.append("-m Project ch0 phiStarll_quantile")

        ### block for variations
        # as
        args.append("--varNames pdfAlphaS --varLabel '$\\alpha_\\mathrm{S}{\\pm}1\\sigma$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # TMP beamfuncs
        # args.append("--varNames resumTNP_b_qg resumTNP_b_qqbarV resumTNP_b_qqDS resumTNP_b_qqS resumTNP_b_qqV --varLabel 'TNP $B_{qg}$' 'TNP $B_{qqV}$' 'TNP $B_{q\\bar{q}V}$' 'TNP $B_{qqS}$' 'TNP $B_{qq\\Delta S}$' --lowerLegCols 3 --lowerLegPos 'upper right'")

        # TNP others
        # args.append("--varNames resumTNP_gamma_cusp resumTNP_gamma_mu_q resumTNP_gamma_nu resumTNP_h_qqV resumTNP_s --varLabel 'TNP $\\Gamma_{\\mathrm{cusp}}$' 'TNP $\\gamma_{\\mu}$' 'TNP $\\gamma_{\\nu}$' 'TNP $H_{qqV}$' 'TNP $S$' --lowerLegCols 3 --lowerLegPos 'upper right'")

        # resum transition
        # args.append("--varNames resumTransitionZSymAvg resumTransitionZSymDiff resumFOScaleZSymAvg resumFOScaleZSymDiff --varLabel 'Match (avg.)' 'Match (diff.)' '$\\mu_{r/f}$ (avg.)' '$\\mu_{r/f}$ (diff.)' --lowerLegCols 3 --lowerLegPos 'upper left'")

        # resum transition
        # args.append("--varNames QCDscaleZinclusive_PtV0_13000helicity_0_SymAvg QCDscaleZinclusive_PtV0_13000helicity_2_SymAvg --varLabel 'P0 (inc)' 'P2 (inc)' --lowerLegCols 3 --lowerLegPos 'upper left'" )

        # New NP: gamma (eigvars)
        # args.append("--varNames scetlibNPgammaEigvar1 scetlibNPgammaEigvar2 scetlibNPgammaEigvar3 --varLabel '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 1' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 2' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 3' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # New NP: gamma (inflated x10)
        # args.append("--varNames scetlibNPgammaLambda2 scetlibNPgammaLambda4 scetlibNPgammaLambdaInf --varLabel '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{2}$' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{4}$' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{\\infty}$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # Old NP: single gamma var
        # args.append("--varNames scetlibNPgamma --varLabel 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # New NP: f
        # args.append("--varNames chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 chargeVgenNP0scetlibNPZdelta_lambda2 --varLabel '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{2}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{4}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Delta\\Lambda^{2}$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # Old NP: f
        # args.append("--varNames chargeVgenNP0scetlibNPZLambda2 chargeVgenNP0scetlibNPZLambda4 chargeVgenNP0scetlibNPZDelta_Lambda2 --varLabel '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{2}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{4}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Delta\\Lambda^{2}$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # PDF: top 3 CT18Z
        # args.append("--varNames pdf14CT18ZSymAvg pdf26CT18ZSymAvg pdf27CT18ZSymAvg --varLabel 'PDF (14)' 'PDF (26)' 'PDF (27)' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # PDF: top 3 MSHT20
        # args.append("--varNames pdf13MSHT20SymAvg pdf16MSHT20SymAvg pdf31MSHT20SymAvg --varLabel 'PDF (13)' 'PDF (16)' 'PDF (31)' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # CT18Z (2) + aS
        # args.append("--varNames pdf26CT18ZSymAvg pdf14CT18ZSymAvg pdfAlphaS --varLabel 'PDF (26)' 'PDF (14)' '$\\alpha_\\mathrm{S}$' --lowerLegCols 2 --lowerLegPos 'upper right'")

        # # overlay lattice prediction (hardcoded external file)
        # arg_variations = " --addVarFromExternalFile /scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints/fitresults_asimov.hdf5 red 'new NP model' --externalRatio num --lowerLegCols 2 --lowerLegPos 'upper right'"


        if prefit:
            args.append("--prefit")
        
        if show_data:
            if not prefit:
                raise ValueError("Trying to plot unblinded postfit distribution, are you sure you want to?")
            args.append("--dataHist data_obs --chisq none")
        else:
            args.append("--noData")

        args = [a for a in args if a]
        all_args = " ".join(args)
        
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_hists.py {all_args}"
        print(command)
        os.system(command)

   
    def plot_impacts(self):


        args = [f"'{os.path.join(self.outfolder_fit, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}'"]
        args.append(f"--title CMS")
        args.append(f"--subtitle WiP")
        args.append(f"--postfix {self.postfix_fit_plot}")
        args.append(f"--showNumbers")
        args.append(f"--diffPullAsym")
        args.append(f"--pullrange '2.1'")
        args.append(f"--config 'utilities/styles/styles.py'")
        args.append(f"--oneSidedImpacts")
        args.append(f"--grouping min")
        args.append(f"--otherExtensions pdf png")
        args.append(f"-n 50")
        args.append(f"--scaleImpacts 2.0")
        args.append(f"-o '{os.path.join(self.outfolder_plots, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'")

        args = [a for a in args if a]
        all_args = " ".join(args)

        # traditional impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py {all_args}"
        print(command)
        os.system(command)
        
        # global impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py {all_args} --globalImpacts"
        print(command)
        os.system(command)



if __name__ == "__main__":
    main()





