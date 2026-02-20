import os
from dataclasses import dataclass, field
from enum import Enum









def main():


    outfolder_base = "/scratch/submit/cms/areimers/wmass"
    plotfolder_base = "/work/submit/areimers/wmass/plots"
    wremnants_folder = "/work/submit/areimers/wmass/WRemnants"

    max_threads = 300

    corrs_and_postfixes_and_legendnames = [
        # ("scetlib_dyturbo", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO\n(nominal)"),
        # ("scetlib_dyturbo", "FixedAlphaS_AsAndPdfFromHels_Inputs4d", "N$^{3{+}0}$LL+NNLO\n(nominal)"),

        # Theory corrections with new NP model. 
        # 1) that exactly reproduces the old central values and approximately the old variations. 
        # 2) With lattice-motivated central values and lattice-constrained uncertainties
        # 3) With lattice-motivated central values and old-like uncertainties

        # ("scetlib_dyturbo_NewNPModel_OldValsAndVars", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, old vals/vars)"),
        # ("scetlib_dyturbo_NewNPModel_LatticeValsOldVars", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/old vars)"),
        # ("scetlib_dyturbo_NewNPModel_LatticeValsAndVars", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars)"),

        # ("scetlib_dyturbo_NewNPModel_LatticeValsOldVars_PdfAsRerun", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/old vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"),
        # ("scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"),
        # ("scetlib_dyturboN3p0LL_LatticeNP", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"), # a check for the committed version
        # ("scetlib_dyturboN3p0LL_LatticeNP", "FixedAlphaS_HackedPDFVarsFromHels", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"), # a check for the committed version

        # ("scetlib_dyturboN3p0LL_LatticeNP_RelToMinnlo", "FixedAlphaS", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"), #

        # ("scetlib_dyturboN3p0LL_LatticeNP", "FixedAlphaS_AsAndPdfFromHels_Inputs4d", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"), 

        ("scetlib_dyturboN3p0LL_LatticeNP", "FixedAlphaS_AsAndPdfFromHels_Inputs4d", "N$^{3{+}0}$LL+NNLO (new NP, lattice vals/vars, PDF+$\\alpha_\\mathrm{s}$ rerun)"), 

    ]



    # modes = [Mode.as_dilepton]
    modes = [Mode.as_dilepton_4d]

    use_own_pdfvars = True
    # use_own_pdfvars = False

    add_oldnp_pdf_as_vars = True
    # add_oldnp_pdf_as_vars = False

    # use_oldnp_pdf_vars = True
    use_oldnp_pdf_vars = False

    # use_oldnp_as_vars = True
    use_oldnp_as_vars = False




    for mode in modes:
        for (theory_corr, postfix, legname) in corrs_and_postfixes_and_legendnames:
            runner = AnalysisRunner(mode=mode, max_threads=max_threads, outfolder_base=outfolder_base, plotfolder_base=plotfolder_base, wremnants_folder=wremnants_folder, theory_corr=theory_corr, use_own_pdfvars=use_own_pdfvars, add_oldnp_pdf_as_vars=add_oldnp_pdf_as_vars, use_oldnp_pdf_vars=use_oldnp_pdf_vars, use_oldnp_as_vars=use_oldnp_as_vars, postfix=postfix)

            # runner.run_histmaker()

            # runner.run_setup_rabbit()
            # runner.run_setup_rabbit(np_model="lattice_oldvars_singlegammavars")
            # runner.run_setup_rabbit(theory_corr_from_hels=True)
            # runner.run_setup_rabbit(theory_corr_from_hels=True, np_model="LatticeNoConstraints")

            # runner.run_fit_rabbit()
            # runner.run_fit_rabbit(np_model="lattice_oldvars_singlegammavars")
            # runner.run_fit_rabbit(theory_corr_from_hels=True)
            # runner.run_fit_rabbit(theory_corr_from_hels=True, np_model="LatticeNoConstraints")

            # runner.plot_fit_rabbit(prefit=True)
            # runner.plot_fit_rabbit(prefit=False) # postfit
            # runner.plot_fit_rabbit(prefit=True, theory_corr_from_hels=True)
            # runner.plot_fit_rabbit(prefit=False, theory_corr_from_hels=True) # postfit
            # runner.plot_fit_rabbit(prefit=True, theory_corr_from_hels=True, np_model="LatticeNoConstraints")
            # runner.plot_fit_rabbit(prefit=False, theory_corr_from_hels=True, np_model="LatticeNoConstraints") # postfit

            # runner.plot_impacts()
            # runner.plot_impacts(np_model="lattice_oldvars_singlegammavars")
            # runner.plot_impacts(theory_corr_from_hels=True)
            # runner.plot_impacts(theory_corr_from_hels=True, np_model="LatticeNoConstraints")

            # runner.plot_param_corr(theory_corr_from_hels=True)
            # runner.plot_param_corr(theory_corr_from_hels=True, np_model="LatticeNoConstraints")





class Mode(Enum):
    as_dilepton = "as_dilepton"
    as_dilepton_4d = "as_dilepton_4d"

mode_configs = {
    Mode.as_dilepton: {
        "exec_histmaker": "mz_dilepton.py",
        "args_histmaker": " --axes yll ptll --csVarsHist",
        "outputfile_histmaker": "mz_dilepton",
        "subfolder_histmaker": "ForAlphaS/WRemDev",
        "args_setup_rabbit": " --fitvar ptll-yll --noi alphaS",
        "fitinput_folder": "ZMassDilepton_ptll_yll",
        "fitinput_filename": "ZMassDilepton.hdf5",
        "args_unblind": "",
    },
}
mode_configs[Mode.as_dilepton_4d] = mode_configs[Mode.as_dilepton].copy()
mode_configs[Mode.as_dilepton_4d]["args_setup_rabbit"] = " --fitvar ptll-yll-cosThetaStarll_quantile-phiStarll_quantile --noi alphaS"
mode_configs[Mode.as_dilepton_4d]["fitinput_folder"] = "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile"


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
    theory_corr: str
    use_own_pdfvars: bool
    add_oldnp_pdf_as_vars: bool
    use_oldnp_pdf_vars: bool
    use_oldnp_as_vars: bool
    postfix: str = ""
    theory_corrs: str = ""
    max_threads: int = -1
    np_model: str = "Delta_Lambda"
    
    arg_pseudodata: str = field(init=False)
    nominal_corr: str = field(init=False)
    postfix_histmaker: str = field(init=False)
    postfix_fit_plot: str = field(init=False)
    exec_histmaker: str = field(init=False)
    args_histmaker: str = field(init=False)
    outputfile_histmaker: str = field(init=False)
    args_setup_rabbit: str = field(init=False)
    fitinput_folder: str = field(init=False)
    fitinput_filename: str = field(init=False)
    args_unblind: str = field(init=False)
    pdfvars_from_corr: str = field(init=False)
    asvars_from_corr: str = field(init=False)

    def __post_init__(self):
        self.theory_corrs = f"'{self.theory_corr}' '{self.theory_corr}_pdfas'" if self.theory_corr else "'scetlib_dyturbo' 'scetlib_dyturboCT18Z_pdfas'"
        self.nominal_corr = self.theory_corrs.split(" ")[0].strip("'")
        if self.use_own_pdfvars:
            self.theory_corrs += f" '{self.theory_corr}_CT18ZVars'"
        if self.theory_corr == "scetlib_dyturbo":
            self.theory_corrs = self.theory_corrs.replace("scetlib_dyturbo_pdfas", "scetlib_dyturboCT18Z_pdfas")
            if self.use_own_pdfvars:
                self.theory_corrs = self.theory_corrs.replace("scetlib_dyturbo_CT18ZVars", "scetlib_dyturboCT18ZVars")
        if "scetlib_dyturbo" not in self.theory_corrs:
            self.theory_corrs += " 'scetlib_dyturbo' 'scetlib_dyturboCT18Z_pdfas'"
        if self.add_oldnp_pdf_as_vars:
            self.theory_corrs += " 'scetlib_dyturboCT18Z_pdfas' 'scetlib_dyturboCT18ZVars'"

        if self.use_own_pdfvars:
            if self.use_oldnp_pdf_vars:
                self.pdfvars_from_corr = "scetlib_dyturboCT18ZVarsCorr"
            else:
                self.pdfvars_from_corr = f"{self.theory_corr}_CT18ZVarsCorr"
        else:
            self.pdfvars_from_corr = ""

        if self.use_oldnp_as_vars:
            self.asvars_from_corr = "scetlib_dyturboCT18Z_pdfasCorr"
        else:
            self.asvars_from_corr = f"{self.theory_corr}_pdfasCorr"

        if self.use_own_pdfvars:
            self.postfix += "_OwnPDFVars"
        if self.add_oldnp_pdf_as_vars:
            self.postfix += "_InclOldNPPDFAsVars"
        
        self.postfix_histmaker = f"{self.postfix}"
        self.postfix_fit_plot = f"{self.nominal_corr}Corr"
        if self.postfix:
            self.postfix_fit_plot += f"_{self.postfix}"
        self.arg_pseudodata = f" --pseudoData '{self.nominal_corr}Corr'" if "scetlib_dyturbo" not in self.nominal_corr else " --pseudoData 'uncorr'"

        if "NewNPModel" in self.theory_corr:
            self.np_model = "tanh2"

            if "LatticeValsAndVars" in self.theory_corr:
                self.np_model = "lattice"

            if "LatticeValsOldVars" in self.theory_corr:
                self.np_model = "lattice_oldvars"

        if "_LatticeNP" in self.theory_corr:
                self.np_model = "LatticeEigvars"
        


        self.exec_histmaker = mode_configs[self.mode]["exec_histmaker"]
        self.args_histmaker = mode_configs[self.mode]["args_histmaker"]
        self.outputfile_histmaker = mode_configs[self.mode]["outputfile_histmaker"]
        self.subfolder_histmaker = mode_configs[self.mode]["subfolder_histmaker"]
        self.args_setup_rabbit = mode_configs[self.mode]["args_setup_rabbit"]
        self.fitinput_folder = mode_configs[self.mode]["fitinput_folder"]
        self.fitinput_filename = mode_configs[self.mode]["fitinput_filename"]
        self.args_unblind = mode_configs[self.mode]["args_unblind"]
        







        


    def run_histmaker(self):

        arg_maxthreads = f" -j {self.max_threads}" if self.max_threads > 0 else ""

        # produces {outfolder_base}/histmaker/mz_dilepton_{theory_corr}Corr_maxFiles_m1_{postfix}.hdf5
        if self.subfolder_histmaker:
            outfolder = os.path.join(self.outfolder_base, 'histmaker', self.subfolder_histmaker)
        else:
            outfolder = os.path.join(self.outfolder_base, 'histmaker')

            
        command = f"python scripts/histmakers/{self.exec_histmaker} -o '{outfolder}' --theoryCorr {self.theory_corrs} {self.args_histmaker} --maxFiles '-1' -p {self.postfix_histmaker}{arg_maxthreads}"
        # command = f"python scripts/histmakers/{self.exec_histmaker} -o '{outfolder}' --theoryCorr {self.theory_corrs} {self.args_histmaker} --maxFiles '1' -p {self.postfix_histmaker}{arg_maxthreads}"

        print(command)
        os.system(command)



    def run_setup_rabbit(self, np_model=None, theory_corr_from_hels=False):
        # produces {outfolder_base}/fitinputs/ZMassDilepton_ptll_{theoryCorr}Corr_{arg_postfix}.hdf5

        arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
        if theory_corr_from_hels:
            arg_theory_corr_hels = ""
            arg_postfix += "_TheoryCorrsFromHels"
        else:
            arg_theory_corr_hels = " --noTheoryCorrsViaHelicities"
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'histmaker', self.subfolder_histmaker)
            outfolder = os.path.join(self.outfolder_base, 'fitinputs', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'histmaker')
            outfolder = os.path.join(self.outfolder_base, 'fitinputs')

        if np_model:
            arg_np_model = np_model
            if np_model == "LatticeNoConstraints":
                arg_postfix += f"_{np_model}"
            else:
                arg_postfix += f"_NPModel{np_model}"
        else:
            arg_np_model = self.np_model

        setup_rabbit_args = self.args_setup_rabbit
        if self.use_own_pdfvars and self.pdfvars_from_corr:
            setup_rabbit_args += f" --pdfUncFromCorr {self.pdfvars_from_corr}"

        if self.asvars_from_corr:
            setup_rabbit_args += f" --asUncFromCorr {self.asvars_from_corr}"

        if self.use_oldnp_pdf_vars:
            arg_postfix += f"_PDFVarsFromOldNP"

        if self.use_oldnp_as_vars:
            arg_postfix += f"_ASVarsFromOldNP"


        command = f"python scripts/rabbit/setupRabbit.py --verbose 3 -i '{os.path.join(infolder, f'{self.outputfile_histmaker}_{self.theory_corr}Corr_maxFiles_m1_{self.postfix_histmaker}.hdf5')}' --npUnc '{arg_np_model}' -o '{outfolder}'{arg_theory_corr_hels} --realData{setup_rabbit_args} -p {arg_postfix}{self.arg_pseudodata}"

        print(command)
        os.system(command)

    def run_fit_rabbit(self, np_model=None, theory_corr_from_hels=False):
        # produces {outfolder_base}/fitresults/ZMassDilepton_ptll_{arg_postfix}/fitresults_data.hdf5

        arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
        if theory_corr_from_hels:
            arg_postfix += "_TheoryCorrsFromHels"

        if np_model:
            if np_model == "LatticeNoConstraints":
                arg_postfix += f"_{np_model}"
            else:
                arg_postfix += f"_NPModel{np_model}"

        if self.use_oldnp_pdf_vars:
            arg_postfix += f"_PDFVarsFromOldNP"
        if self.use_oldnp_as_vars:
            arg_postfix += f"_ASVarsFromOldNP"
            
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitinputs', self.subfolder_histmaker)
            outfolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitinputs')
            outfolder = os.path.join(self.outfolder_base, 'fitresults')


        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', self.fitinput_filename)}' --globalImpacts --computeHistErrors --computeVariations --saveHists --saveHistsPerProcess --computeHistCov -m Basemodel -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile --doImpacts -t -1 --postfix asimov -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}'{self.args_unblind}"
        print(command)
        os.system(command)


        # command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', self.fitinput_filename)}' --globalImpacts --computeHistErrors --computeVariations --saveHists --saveHistsPerProcess --computeHistCov -m Basemodel -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile --doImpacts -t 0 --postfix data -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}'{self.args_unblind}"
        # print(command)
        # os.system(command)




    def plot_fit_rabbit(self, np_model=None, theory_corr_from_hels=False, prefit=True):
        # produces {outfolder_base}/fitresults/ZMassDilepton_ptll_{arg_postfix}/fitresults_data.hdf5

        arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
        if theory_corr_from_hels:
            arg_postfix += "_TheoryCorrsFromHels"

        if np_model:
            if np_model == "LatticeNoConstraints":
                arg_postfix += f"_{np_model}"
            else:
                arg_postfix += f"_NPModel{np_model}"

        if self.use_oldnp_pdf_vars:
            arg_postfix += f"_PDFVarsFromOldNP"
        if self.use_oldnp_as_vars:
            arg_postfix += f"_ASVarsFromOldNP"
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
            infolder_histmaker = os.path.join(self.outfolder_base, 'histmaker', self.subfolder_histmaker)
            outfolder = os.path.join(self.plotfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitresults')
            infolder_histmaker = os.path.join(self.outfolder_base, 'histmaker')
            outfolder = os.path.join(self.plotfolder_base, 'fitresults')

        arg_prefit = " --prefit" if prefit else ""

        if prefit:
            arg_rrange = " --rrange 0.85 1.15"
        else: 
            arg_rrange = " --rrange 0.95 1.05"


        if self.mode in [Mode.as_dilepton, Mode.as_dilepton_4d]:
            # no var
            # arg_variations = ""
            
            # alpha_s var
            # arg_variations = " --varNames pdfAlphaS --varLabel '$\\alpha_\\mathrm{S}{\\pm}1\\sigma$' --lowerLegCols 2 --lowerLegPos 'upper right'"

            # # overlay lattice prediction (hardcoded external file)
            # arg_variations = " --addVarFromExternalFile /scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_FixedAlphaS_AsAndPdfFromHels_Inputs4d_AlphaS_TheoryCorrsFromHels/fitresults_asimov.hdf5 red 'new NP model' --externalRatio num --lowerLegCols 2 --lowerLegPos 'upper right'"

            # strongest TNP var
            # arg_variations = " --varNames resumTNP_gamma_nu --varLabel 'TNP $\\gamma_\\nu$'"

            # Gamma NP vars
            if self.theory_corr in ["scetlib_dyturbo_NewNPModel_LatticeValsAndVars", "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun", "scetlib_dyturboN3p0LL_LatticeNP"]:
                if np_model == "LatticeNoConstraints":
                    # arg_variations = " --varNames scetlibNPgammalambda2nu scetlibNPgammalambda4nu scetlibNPgammalambdainfnu --varLabel 'SCETLib $\\gamma$: $\\lambda^{2}_{\\nu}$' 'SCETLib $\\gamma$: $\\lambda^{4}_{\\nu}$' 'SCETLib $\\gamma$: $\\lambda^{\\infty}_{\\nu}$' --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
                    arg_variations = " --varNames scetlibNPgammaLambda2 scetlibNPgammaLambda4 scetlibNPgammaLambdaInf --varLabel 'NP $\\gamma$: $\\lambda^{2}_{\\nu}$' 'NP $\\gamma$: $\\lambda^{4}_{\\nu}$' 'NP: $\\gamma$: $\\lambda^{\\infty}_{\\nu}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"
                else:
                    arg_variations = " --varNames scetlibNPgammaEigvar1 scetlibNPgammaEigvar2 scetlibNPgammaEigvar3 --varLabel 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$ eig. var. 1' 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$ eig. var. 2' 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$ eig. var. 3' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"
            else:
                arg_variations = " --varNames scetlibNPgamma --varLabel 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

            # # Lambda NP vars
            # if "NewNPModel" in self.theory_corr or "_LatticeNP" in self.theory_corr:
            #     arg_variations = " --varNames chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 chargeVgenNP0scetlibNPZdelta_lambda2 --varLabel 'SCETLib $\\Lambda^{2}$' 'SCETLib $\\Lambda^{4}$' 'SCETLib $\\Delta\\Lambda^{2}$' --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
            # else:
            #     arg_variations = " --varNames chargeVgenNP0scetlibNPZLambda2 chargeVgenNP0scetlibNPZLambda4 chargeVgenNP0scetlibNPZDelta_Lambda2 --varLabel 'SCETLib $\\Lambda^{2}$' 'SCETLib $\\Lambda^{4}$' 'SCETLib $\\Delta\\Lambda^{2}$' --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"

            # some of the largest PDF vars
            # arg_variations = " --varNames pdf14CT18ZSymAvg pdf26CT18ZSymAvg pdf27CT18ZSymAvg --varLabel 'PDF (14)' 'PDF (26)' 'PDF (27)' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"
            # arg_variations = " --varNames pdf16CT18ZSymAvg pdf18CT18ZSymAvg pdf25CT18ZSymAvg --varLabel 'PDF (16)' 'PDF (18)' 'PDF (25)' --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
            # allpdfs_nuisnames = " ".join([f"pdf{i}CT18ZSymDiff" for i in range(1,30)])
            # allpdfs_nicenames = " ".join([f"'PDF ({i})'" for i in range(1,30)])
            # arg_variations = f" --varNames {allpdfs_nuisnames} --varLabel {allpdfs_nicenames} --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
            if prefit:
                arg_rrange = " --rrange 0.95 1.05"
                # arg_rrange = " --rrange 0.97 1.03"
                # arg_rrange = " --rrange 0.85 1.15"

            # some of the largest TNP vars
            # arg_variations = " --varNames resumTNP_b_qqV resumTNP_s resumTNP_b_qg --varLabel 'TNP BF qqV' 'TNP soft func.' 'TNP BF qg' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

            # PDFs with largest change in global impact + pdfas variation
            # arg_variations = " --varNames pdf26CT18ZSymAvg pdf14CT18ZSymAvg pdfAlphaS --varLabel 'PDF (26)' 'PDF (14)' '$\\alpha_\\mathrm{S}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

        arg_projections = " -m Project ch0 ptll"
        # arg_projections = " -m Project ch0 ptll -m Project ch0 yll"
        # arg_projections = " -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile"
        #  --noData
        # --dataHist data_obs --chisq none

        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_hists.py --config utilities/styles/styles.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}'  --title CMS --titlePos 1 --subtitle WiP --extraTextLoc 0.5 0.94 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}' --legCols 1 --yscale 1.2{arg_prefit} --legCol 1 --noData --legSize small{arg_projections}{arg_variations}{arg_rrange}"
        print(command)
        os.system(command)

        # if prefit:
        #     command = f"python scripts/plotting/makeDataMCStackPlot.py --ratioToData --yscale 1.2 --baseName nominal --hists ptll -o /work/submit/areimers/wmass/plots -f '{os.path.join('distributions', f'{self.fitinput_folder}_{arg_postfix}')}' -p {arg_postfix} '{os.path.join(infolder_histmaker, f'{self.outputfile_histmaker}_{self.theory_corr}Corr_maxFiles_m1_{self.postfix_histmaker}.hdf5')}' variation --varName nominal_uncorr --varLabel MiNNLO --colors red"
        #     print(command)
        #     os.system(command)

   
    def plot_impacts(self, np_model=None, theory_corr_from_hels=False):

        arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
        if theory_corr_from_hels:
            arg_postfix += "_TheoryCorrsFromHels"

        if np_model:
            if np_model == "LatticeNoConstraints":
                arg_postfix += f"_{np_model}"
            else:
                arg_postfix += f"_NPModel{np_model}"

        if self.use_oldnp_pdf_vars:
            arg_postfix += f"_PDFVarsFromOldNP"
        if self.use_oldnp_as_vars:
            arg_postfix += f"_ASVarsFromOldNP"
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
            outfolder = os.path.join(self.plotfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitresults')
            outfolder = os.path.join(self.plotfolder_base, 'fitresults')

        # traditional impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}' --title CMS --subtitle WiP --postfix {arg_postfix} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping min --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}'"
        print(command)
        os.system(command)

        # global impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}' --globalImpacts --title CMS --subtitle WiP --postfix {arg_postfix} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping min --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}'"
        print(command)
        os.system(command)

        # pulls
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}' --title CMS --subtitle WiP --postfix {arg_postfix} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping max --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}' --noImpacts"
        print(command)
        os.system(command)
        

    def plot_param_corr(self, np_model=None, theory_corr_from_hels=False):

        arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
        if theory_corr_from_hels:
            arg_postfix += "_TheoryCorrsFromHels"

        if np_model:
            if np_model == "LatticeNoConstraints":
                arg_postfix += f"_{np_model}"
            else:
                arg_postfix += f"_NPModel{np_model}"

        if self.use_oldnp_pdf_vars:
            arg_postfix += f"_PDFVarsFromOldNP"
        if self.use_oldnp_as_vars:
            arg_postfix += f"_ASVarsFromOldNP"
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
            outfolder = os.path.join(self.plotfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitresults')
            outfolder = os.path.join(self.plotfolder_base, 'fitresults')

        if self.theory_corr in ["scetlib_dyturbo_NewNPModel_LatticeValsAndVars", "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun", "scetlib_dyturboN3p0LL_LatticeNP"]:
            arg_nplambda = "chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 chargeVgenNP0scetlibNPZdelta_lambda2"
            if np_model == "LatticeNoConstraints":
                arg_npgamma = "scetlibNPgammaLambda2 scetlibNPgammaLambda4 scetlibNPgammaLambdaInf"
            else:
                arg_npgamma = "scetlibNPgammaEigvar1 scetlibNPgammaEigvar2 scetlibNPgammaEigvar3"
        else:
            arg_npgamma = "scetlibNPgamma"
            arg_nplambda = "chargeVgenNP0scetlibNPZLambda2 chargeVgenNP0scetlibNPZLambda4 chargeVgenNP0scetlibNPZDelta_Lambda2"

        command = f"/work/submit/areimers/wmass/WRemnants/rabbit/bin/rabbit_plot_hists_cov.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}' --result asimov --config 'utilities/styles/styles.py' --correlation --showNumbers --title CMS --subtitle WiP --titlePos 0 --scaleTextSize '0.4' -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}' --params {arg_npgamma} {arg_nplambda} 'resumTNP_gamma_nu' lumi pdf14CT18ZSymAvg pdf26CT18ZSymAvg"
        print(command)
        os.system(command)


if __name__ == "__main__":
    main()





