import os
from dataclasses import dataclass, field
from enum import Enum









def main():


    outfolder_base = "/scratch/submit/cms/areimers/wmass"
    plotfolder_base = "/work/submit/areimers/wmass/plots"
    wremnants_folder = "/work/submit/areimers/wmass/WRemnants"

    max_threads = 350

    configs = [

        ##### Old NP
        # {# CT18Z: PDF vars from sc+dy
        #     "theory_corrs": ["scetlib_dyturbo", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "OldNP_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturboCT18ZVarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboCT18Z_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "OldNP_AllCT18Z_PDFFromCorr_ASFromCorr",
        #     "np_model": "Delta_Lambda"
        # },
        # {# MSHT20: PDF vars from sc+dy
        #     "theory_corrs": ["scetlib_dyturboMSHT20", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "OldNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboMSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboMSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "OldNP_AllMSHT20_PDFFromCorr_ASFromCorr",
        #     "np_model": "Delta_Lambda"
        # },
        # {# MSHT20: PDF vars from sc+dy, inflated x1.70
        #     "theory_corrs": ["scetlib_dyturboMSHT20", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "OldNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboMSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboMSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.70,
        #     "postfix_rabbit": "OldNP_AllMSHT20_PDFFromCorr_ASFromCorr_PDFInflated1p7",
        #     "np_model": "Delta_Lambda"
        # },

        # {# CT18Z: PDF vars from minnlo
        #     "theory_corrs": ["scetlib_dyturbo", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "OldNP_AllCT18Z",
        #     "pdfvars_from_corr": "",
        #     "asvars_from_corr": "scetlib_dyturboCT18Z_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "OldNP_AllCT18Z_PDFFromMinnlo_ASFromCorr",
        #     "np_model": "Delta_Lambda"
        # },
        # {# MSHT20: PDF vars from minnlo
        #     "theory_corrs": ["scetlib_dyturboMSHT20", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "OldNP_AllMSHT20",
        #     "pdfvars_from_corr": "",
        #     "asvars_from_corr": "scetlib_dyturboMSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "OldNP_AllMSHT20_PDFFromMinnlo_ASFromCorr",
        #     "np_model": "Delta_Lambda"
        # },
        # {# MSHT20: PDF vars from minnlo
        #     "theory_corrs": ["scetlib_dyturboMSHT20", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "OldNP_AllMSHT20",
        #     "pdfvars_from_corr": "",
        #     "asvars_from_corr": "scetlib_dyturboMSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.70,
        #     "postfix_rabbit": "OldNP_AllMSHT20_PDFFromMinnlo_ASFromCorr_PDFInflated1p7",
        #     "np_model": "Delta_Lambda"
        # },

        
        ##### New NP with lattice-constrained eigenvariations
        {
            "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP", "scetlib_dyturboN3p0LL_LatticeNP_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVars", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
            "pdfs": ["ct18z"],
            "postfix_histmaker": "NewNP_AllCT18Z",
            "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr",
            "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr",
            "pdf_inflation_factor": 1.0,
            "postfix_rabbit": "NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints",
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
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP_MSHT20", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20Vars", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "NewNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_AllMSHT20_PDFFromCorr_ASFromCorr_WithLatticeConstraints",
        #     "np_model": "LatticeEigvars"
        # },


        # ##### New NP, no lattice constraints, instead inflating the individual vars x20
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP", "scetlib_dyturboN3p0LL_LatticeNP_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVars", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_CT18ZVarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_NoLatticeConstraints",
        #     "np_model": "LatticeNoConstraints"
        # },
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP_MSHT20", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20Vars", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "NewNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_AllMSHT20_PDFFromCorr_ASFromCorr_NoLatticeConstraints",
        #     "np_model": "LatticeNoConstraints"
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
        # {
        #     "theory_corrs": ["scetlib_dyturboN3p0LL_LatticeNP_MSHT20", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfas", "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20Vars", "scetlib_dyturboMSHT20_pdfas", "scetlib_dyturboMSHT20Vars"],
        #     "pdfs": ["msht20"],
        #     "postfix_histmaker": "NewNP_AllMSHT20",
        #     "pdfvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_MSHT20VarsCorr",
        #     "asvars_from_corr": "scetlib_dyturboN3p0LL_LatticeNP_MSHT20_pdfasCorr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_AllMSHT20_PDFFromCorr_ASFromCorr_NoLatticeConstraintsx10",
        #     "np_model": "LatticeNoConstraintsx10"
        # },

        
        ##### New NP in FINE BINNING and two bins in Q=[10, 60, 120] with the two lattice models (constrained and x10)
        # {
        #     "theory_corrs": ["scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_", "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfas_", "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfvars_", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_FineBins_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfvars_Corr",
        #     "asvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfas_Corr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_FineBins_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints",
        #     "np_model": "LatticeEigvars"
        # },
        # {
        #     "theory_corrs": ["scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_", "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfas_", "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfvars_", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_FineBins_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfvars_Corr",
        #     "asvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_CT18Z_N3p0LL_N2L0_pdfas_Corr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_FineBins_AllCT18Z_PDFFromCorr_ASFromCorr_NoLatticeConstraintsx10",
        #     "np_model": "LatticeNoConstraintsx10"
        # },

        
        ##### New NP in FINE BINNING and ONE bin in Q=[60, 120] (hacked make_theory_corr) with the two lattice models (constrained and x10)
        # {
        #     "theory_corrs": ["scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_", "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfas_", "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfvars_", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_FineBins_Q60to120_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfvars_Corr",
        #     "asvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfas_Corr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_FineBins_Q60to120_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints",
        #     "np_model": "LatticeEigvars"
        # },
        # {
        #     "theory_corrs": ["scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_", "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfas_", "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfvars_", "scetlib_dyturboCT18Z_pdfas", "scetlib_dyturboCT18ZVars"],
        #     "pdfs": ["ct18z"],
        #     "postfix_histmaker": "NewNP_FineBins_Q60to120_AllCT18Z",
        #     "pdfvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfvars_Corr",
        #     "asvars_from_corr": "scetlib_dyturbo_LatticeNP_FineBins_Q60to120_CT18Z_N3p0LL_N2L0_pdfas_Corr",
        #     "pdf_inflation_factor": 1.0,
        #     "postfix_rabbit": "NewNP_FineBins_Q60to120_AllCT18Z_PDFFromCorr_ASFromCorr_NoLatticeConstraintsx10",
        #     "np_model": "LatticeNoConstraintsx10"
        # },
    ]



    # modes = [Mode.as_dilepton]
    modes = [Mode.as_dilepton_4d]
    # modes = [Mode.as_dilepton_1d]




    for mode in modes:
        for config in configs:
            runner = AnalysisRunner(mode=mode, max_threads=max_threads, outfolder_base=outfolder_base, plotfolder_base=plotfolder_base, wremnants_folder=wremnants_folder, theory_corrs=config["theory_corrs"], pdfs=config["pdfs"], pdfvars_from_corr=config["pdfvars_from_corr"], asvars_from_corr=config["asvars_from_corr"], pdf_inflation_factor=config["pdf_inflation_factor"], np_model=config["np_model"], postfix_histmaker=config["postfix_histmaker"], postfix_rabbit=config["postfix_rabbit"])

            # runner.run_histmaker()

            # runner.run_setup_rabbit()

            # runner.run_fit_rabbit()

            runner.plot_fit_rabbit(prefit=True)
            # runner.plot_fit_rabbit(prefit=False) # postfit

            # runner.plot_impacts()

            # runner.plot_param_corr(theory_corr_from_hels=True)
            # runner.plot_param_corr(theory_corr_from_hels=True, np_model="LatticeNoConstraints")





class Mode(Enum):
    as_dilepton = "as_dilepton"
    as_dilepton_4d = "as_dilepton_4d"
    as_dilepton_1d = "as_dilepton_1d"

mode_configs = {
    Mode.as_dilepton: {
        "exec_histmaker": "mz_dilepton.py",
        "args_histmaker": " --axes yll ptll --csVarsHist",
        "outputfile_histmaker": "mz_dilepton",
        "subfolder_histmaker": "ForAlphaS/WRemDev/NewSteer",
        "args_setup_rabbit": " --fitvar ptll-yll --noi alphaS",
        "fitinput_folder": "ZMassDilepton_ptll_yll",
        "fitinput_filename": "ZMassDilepton.hdf5",
        "args_unblind": "",
    },
}
mode_configs[Mode.as_dilepton_4d] = mode_configs[Mode.as_dilepton].copy()
mode_configs[Mode.as_dilepton_4d]["args_setup_rabbit"] = " --fitvar ptll-yll-cosThetaStarll_quantile-phiStarll_quantile --noi alphaS"
mode_configs[Mode.as_dilepton_4d]["fitinput_folder"] = "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile"

mode_configs[Mode.as_dilepton_1d] = mode_configs[Mode.as_dilepton].copy()
mode_configs[Mode.as_dilepton_1d]["args_setup_rabbit"] = " --fitvar ptll --noi alphaS"
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
    pdfvars_from_corr: str
    asvars_from_corr: str
    pdf_inflation_factor: float
    np_model: str
    postfix_histmaker: str
    postfix_rabbit: str
    pdfs: list[str] = field(init=["ct18z"])
    max_threads: int = -1
    
    nominal_corr: str = field(init=False)
    postfix_fit_plot: str = field(init=False)
    exec_histmaker: str = field(init=False)
    args_histmaker: str = field(init=False)
    outputfile_histmaker: str = field(init=False)
    args_setup_rabbit: str = field(init=False)
    fitinput_folder: str = field(init=False)
    fitinput_filename: str = field(init=False)
    args_unblind: str = field(init=False)

    def __post_init__(self):
        # self.theory_corrs = f"'{self.theory_corr}' '{self.theory_corr}_pdfas'" if self.theory_corr else "'scetlib_dyturbo' 'scetlib_dyturboCT18Z_pdfas'"
        self.nominal_corr = self.theory_corrs[0]

        self.postfix_histmaker = f"{self.postfix_histmaker}"
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
        







        


    def run_histmaker(self):

        arg_maxthreads = f" -j {self.max_threads}" if self.max_threads > 0 else ""
        arg_theory_corrs = " ".join(self.theory_corrs)

        # produces {outfolder_base}/histmaker/mz_dilepton_{theory_corr}Corr_maxFiles_m1_{postfix}.hdf5
        if self.subfolder_histmaker:
            outfolder = os.path.join(self.outfolder_base, 'histmaker', self.subfolder_histmaker)
        else:
            outfolder = os.path.join(self.outfolder_base, 'histmaker')

        arg_pdfs = "--pdfs " + " ".join(self.pdfs)

            
        command = f"python scripts/histmakers/{self.exec_histmaker} -o '{outfolder}' --theoryCorr {arg_theory_corrs} {arg_pdfs} {self.args_histmaker} --maxFiles '-1' -p {self.postfix_histmaker}{arg_maxthreads}"

        print(command)
        os.system(command)



    def run_setup_rabbit(self):
        # produces {outfolder_base}/fitinputs/ZMassDilepton_ptll_{theoryCorr}Corr_{arg_postfix}.hdf5

        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'histmaker', self.subfolder_histmaker)
            outfolder = os.path.join(self.outfolder_base, 'fitinputs', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'histmaker')
            outfolder = os.path.join(self.outfolder_base, 'fitinputs')

        setup_rabbit_args = self.args_setup_rabbit
        if self.pdfvars_from_corr:
            setup_rabbit_args += f" --pdfUncFromCorr {self.pdfvars_from_corr}"

        if self.asvars_from_corr:
            setup_rabbit_args += f" --asUncFromCorr {self.asvars_from_corr}"

        postfix_histmaker_arg = self.postfix_histmaker
        if self.pdfs[0] != "ct18z":
            postfix_histmaker_arg = f"{self.pdfs[0]}_{self.postfix_histmaker}"



        command = f"python scripts/rabbit/setupRabbit.py --verbose 3 -i '{os.path.join(infolder, f'{self.outputfile_histmaker}_{self.nominal_corr}Corr_maxFiles_m1_{postfix_histmaker_arg}.hdf5')}' --npUnc '{self.np_model}' -o '{outfolder}' --realData{setup_rabbit_args} --scalePdf {self.pdf_inflation_factor} -p {self.postfix_fit_plot}"

        print(command)
        os.system(command)

    def run_fit_rabbit(self):
        # produces {outfolder_base}/fitresults/ZMassDilepton_ptll_{arg_postfix}/fitresults_data.hdf5            
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitinputs', self.subfolder_histmaker)
            outfolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitinputs')
            outfolder = os.path.join(self.outfolder_base, 'fitresults')


        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', self.fitinput_filename)}' --globalImpacts --computeHistErrors --computeVariations --saveHists --saveHistsPerProcess --computeHistCov -m Basemodel -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile --doImpacts -t -1 --postfix asimov -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'{self.args_unblind}"
        print(command)
        os.system(command)

        # command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', self.fitinput_filename)}' --globalImpacts --computeHistErrors --computeVariations --saveHists --saveHistsPerProcess --computeHistCov -m Basemodel -m Project ch0 ptll --doImpacts -t -1 --postfix asimov -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'{self.args_unblind}"
        # print(command)
        # os.system(command)


        # command = f"{self.wremnants_folder}/rabbit/bin/rabbit_fit.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', self.fitinput_filename)}' --globalImpacts --computeHistErrors --computeVariations --saveHists --saveHistsPerProcess --computeHistCov -m Basemodel -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile --doImpacts -t 0 --postfix data -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}'{self.args_unblind}"
        # print(command)
        # os.system(command)




    def plot_fit_rabbit(self, prefit=True):
        # produces {outfolder_base}/fitresults/ZMassDilepton_ptll_{arg_postfix}/fitresults_data.hdf5
        
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
            arg_variations = ""
            
            # alpha_s var
            # arg_variations = " --varNames pdfAlphaS --varLabel '$\\alpha_\\mathrm{S}{\\pm}1\\sigma$' --lowerLegCols 2 --lowerLegPos 'upper right'"

            # # overlay lattice prediction (hardcoded external file)
            # arg_variations = " --addVarFromExternalFile /scratch/submit/cms/areimers/wmass/fitresults/ForAlphaS/WRemDev/NewSteer/ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_scetlib_dyturboN3p0LL_LatticeNPCorr_NewNP_AllCT18Z_PDFFromCorr_ASFromCorr_WithLatticeConstraints/fitresults_asimov.hdf5 red 'new NP model' --externalRatio num --lowerLegCols 2 --lowerLegPos 'upper right'"

            # TNP beamfuncs
            # arg_variations = " --varNames resumTNP_b_qg resumTNP_b_qqbarV resumTNP_b_qqDS resumTNP_b_qqS resumTNP_b_qqV --varLabel 'TNP $B_{qg}$' 'TNP $B_{qqV}$' 'TNP $B_{q\\bar{q}V}$' 'TNP $B_{qqS}$' 'TNP $B_{qq\\Delta S}$' --lowerLegCols 3 --lowerLegPos 'upper right'" 

            # # resum transition
            # arg_variations = " --varNames resumTransitionZSymAvg resumTransitionZSymDiff resumFOScaleZSymAvg resumFOScaleZSymDiff --varLabel 'Match (avg.)' 'Match (diff.)' '$\\mu_{r/f}$ (avg.)' '$\\mu_{r/f}$ (diff.)' --lowerLegCols 3 --lowerLegPos 'upper left' --rrange 0.90 1.15" 

            # # resum transition
            arg_variations = " --varNames QCDscaleZinclusive_PtV0_13000helicity_0_SymAvg QCDscaleZinclusive_PtV0_13000helicity_2_SymAvg --varLabel 'P0 (inc)' 'P2 (inc)' --lowerLegCols 3 --lowerLegPos 'upper left' --rrange 0.98 1.02" 


            

            # TNP others
            # arg_variations = " --varNames resumTNP_gamma_cusp resumTNP_gamma_mu_q resumTNP_gamma_nu resumTNP_h_qqV resumTNP_s --varLabel 'TNP $\\Gamma_{\\mathrm{cusp}}$' 'TNP $\\gamma_{\\mu}$' 'TNP $\\gamma_{\\nu}$' 'TNP $H_{qqV}$' 'TNP $S$' --lowerLegCols 3 --lowerLegPos 'upper right'" 

            # # Gamma NP vars
            # if self.np_model in ["LatticeNoConstraints", "LatticeNoConstraintsx10"]:
            #     # arg_variations = " --varNames scetlibNPgammalambda2nu scetlibNPgammalambda4nu scetlibNPgammalambdainfnu --varLabel 'SCETLib $\\gamma$: $\\lambda^{2}_{\\nu}$' 'SCETLib $\\gamma$: $\\lambda^{4}_{\\nu}$' 'SCETLib $\\gamma$: $\\lambda^{\\infty}_{\\nu}$' --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
            #     arg_variations = " --varNames scetlibNPgammaLambda2 scetlibNPgammaLambda4 scetlibNPgammaLambdaInf --varLabel '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{2}$' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{4}$' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$: $\\lambda_{\\infty}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right' --rrange 0.80 1.20"
            # elif self.np_model == "LatticeEigvars":
            #     arg_variations = " --varNames scetlibNPgammaEigvar1 scetlibNPgammaEigvar2 scetlibNPgammaEigvar3 --varLabel '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 1' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 2' '$\\tilde{\\gamma}_{\\zeta}^{\\mathrm{NP}}$ eig. var. 3' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right' --rrange 0.80 1.20"
            # else:
            #     arg_variations = " --varNames scetlibNPgamma --varLabel 'NP: $\\tilde{\\gamma}^{\\rm{np}}_{\\zeta}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

            # # Lambda NP vars
            # if self.np_model in ["LatticeNoConstraints", "LatticeNoConstraintsx10", "LatticeEigvars"]:
            #     arg_variations = " --varNames chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 chargeVgenNP0scetlibNPZdelta_lambda2 --varLabel '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{2}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{4}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Delta\\Lambda^{2}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right' --rrange 0.80 1.20"
            # else:
            #     arg_variations = " --varNames chargeVgenNP0scetlibNPZLambda2 chargeVgenNP0scetlibNPZLambda4 chargeVgenNP0scetlibNPZDelta_Lambda2 --varLabel '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{2}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Lambda^{4}$' '$\\tilde{f}_i^{\\mathrm{NP}}: \\Delta\\Lambda^{2}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

            # # # some of the largest PDF vars
            # if self.pdfs[0] == "ct18z":
            #     arg_variations = " --varNames pdf14CT18ZSymAvg pdf26CT18ZSymAvg pdf27CT18ZSymAvg --varLabel 'PDF (14)' 'PDF (26)' 'PDF (27)' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"
            #     # arg_variations = " --varNames pdf3CT18ZSymAvg pdf14CT18ZSymAvg --varLabel 'PDF (3)' 'PDF (14)' --showVariations both --lowerLegCols 3 --lowerLegPos 'upper right'"
            # elif self.pdfs[0] == "msht20":
            #     arg_variations = " --varNames pdf13MSHT20SymAvg pdf16MSHT20SymAvg pdf31MSHT20SymAvg --varLabel 'PDF (13)' 'PDF (16)' 'PDF (31)' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"
            # else:
            #     raise ValueError(f"PDF variations to be plotted not defined for PDF {self.pdfs[0]}.")

            # allpdfs_nuisnames = " ".join([f"pdf{i}CT18ZSymDiff" for i in range(1,30)])
            # allpdfs_nicenames = " ".join([f"'PDF ({i})'" for i in range(1,30)])
            # arg_variations = f" --varNames {allpdfs_nuisnames} --varLabel {allpdfs_nicenames} --showVariations both --lowerLegCols 4 --lowerLegPos 'upper right'"
            if prefit:
                arg_rrange = " --rrange 0.98 1.02"
                # arg_rrange = " --rrange 0.95 1.05"
                # arg_rrange = " --rrange 0.97 1.03"
                # arg_rrange = " --rrange 0.85 1.15"

            # some of the largest TNP vars
            # arg_variations = " --varNames resumTNP_b_qqV resumTNP_s resumTNP_b_qg --varLabel 'TNP BF qqV' 'TNP soft func.' 'TNP BF qg' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

            # PDFs with largest change in global impact + pdfas variation
            # arg_variations = " --varNames pdf26CT18ZSymAvg pdf14CT18ZSymAvg pdfAlphaS --varLabel 'PDF (26)' 'PDF (14)' '$\\alpha_\\mathrm{S}$' --showVariations both --lowerLegCols 2 --lowerLegPos 'upper right'"

        # arg_projections = " -m Project ch0 ptll"
        # arg_projections = " -m Project ch0 ptll -m Project ch0 yll"
        arg_projections = " -m Project ch0 ptll -m Project ch0 yll -m Project ch0 cosThetaStarll_quantile -m Project ch0 phiStarll_quantile"
        #  --noData
        # --dataHist data_obs --chisq none

        # with data
        # command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_hists.py --config utilities/styles/styles.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}'  --title CMS --titlePos 1 --subtitle WiP --extraTextLoc 0.5 0.94 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}' --legCols 1 --yscale 1.2{arg_prefit} --legCol 1 --dataHist data_obs --chisq none --legSize small{arg_projections}{arg_variations} --rrange 0.80 1.20"
        # print(command)
        # os.system(command)

        # without data
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_hists.py --config utilities/styles/styles.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}'  --title CMS --titlePos 1 --subtitle WiP --extraTextLoc 0.5 0.94 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}' --showVariations both --legCols 1 --yscale 1.2{arg_prefit} --legCol 1 --noData --legSize small{arg_projections}{arg_rrange}{arg_variations}"
        print(command)
        os.system(command)

        # if prefit:
        #     command = f"python scripts/plotting/makeDataMCStackPlot.py --ratioToData --yscale 1.2 --baseName nominal --hists ptll -o /work/submit/areimers/wmass/plots -f '{os.path.join('distributions', f'{self.fitinput_folder}_{arg_postfix}')}' -p {arg_postfix} '{os.path.join(infolder_histmaker, f'{self.outputfile_histmaker}_{self.theory_corr}Corr_maxFiles_m1_{self.postfix_histmaker}.hdf5')}' variation --varName nominal_uncorr --varLabel MiNNLO --colors red"
        #     print(command)
        #     os.system(command)

   
    def plot_impacts(self):
        
        if self.subfolder_histmaker:
            infolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
            outfolder = os.path.join(self.plotfolder_base, 'fitresults', self.subfolder_histmaker)
        else:
            infolder = os.path.join(self.outfolder_base, 'fitresults')
            outfolder = os.path.join(self.plotfolder_base, 'fitresults')

        # traditional impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}' --title CMS --subtitle WiP --postfix {self.postfix_fit_plot} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping min --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'"
        print(command)
        os.system(command)

        # global impacts
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}' --globalImpacts --title CMS --subtitle WiP --postfix {self.postfix_fit_plot} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping min --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}'"
        print(command)
        os.system(command)

        # pulls
        command = f"{self.wremnants_folder}/rabbit/bin/rabbit_plot_pulls_and_impacts.py '{os.path.join(infolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}', f'fitresults_asimov.hdf5')}' --title CMS --subtitle WiP --postfix {self.postfix_fit_plot} --showNumbers --diffPullAsym --pullrange '2.1' --config 'utilities/styles/styles.py' --oneSidedImpacts --grouping max --otherExtensions pdf png -n 50 --scaleImpacts 2.0 -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{self.postfix_fit_plot}')}' --noImpacts"
        print(command)
        os.system(command)
        

    # def plot_param_corr(self, np_model=None, theory_corr_from_hels=False):

    #     arg_postfix = self.postfix_fit_plot + "_AlphaS"
        
    #     if theory_corr_from_hels:
    #         arg_postfix += "_TheoryCorrsFromHels"

    #     if np_model:
    #         if np_model == "LatticeNoConstraints":
    #             arg_postfix += f"_{np_model}"
    #         else:
    #             arg_postfix += f"_NPModel{np_model}"

    #     if self.use_oldnp_pdf_vars:
    #         arg_postfix += f"_PDFVarsFromOldNP"
    #     if self.use_oldnp_as_vars:
    #         arg_postfix += f"_ASVarsFromOldNP"
        
    #     if self.subfolder_histmaker:
    #         infolder = os.path.join(self.outfolder_base, 'fitresults', self.subfolder_histmaker)
    #         outfolder = os.path.join(self.plotfolder_base, 'fitresults', self.subfolder_histmaker)
    #     else:
    #         infolder = os.path.join(self.outfolder_base, 'fitresults')
    #         outfolder = os.path.join(self.plotfolder_base, 'fitresults')

    #     if self.theory_corr in ["scetlib_dyturbo_NewNPModel_LatticeValsAndVars", "scetlib_dyturbo_NewNPModel_LatticeValsAndVars_PdfAsRerun", "scetlib_dyturboN3p0LL_LatticeNP"]:
    #         arg_nplambda = "chargeVgenNP0scetlibNPZlambda2 chargeVgenNP0scetlibNPZlambda4 chargeVgenNP0scetlibNPZdelta_lambda2"
    #         if np_model == "LatticeNoConstraints":
    #             arg_npgamma = "scetlibNPgammaLambda2 scetlibNPgammaLambda4 scetlibNPgammaLambdaInf"
    #         else:
    #             arg_npgamma = "scetlibNPgammaEigvar1 scetlibNPgammaEigvar2 scetlibNPgammaEigvar3"
    #     else:
    #         arg_npgamma = "scetlibNPgamma"
    #         arg_nplambda = "chargeVgenNP0scetlibNPZLambda2 chargeVgenNP0scetlibNPZLambda4 chargeVgenNP0scetlibNPZDelta_Lambda2"

    #     command = f"/work/submit/areimers/wmass/WRemnants/rabbit/bin/rabbit_plot_hists_cov.py '{os.path.join(infolder, f'{self.fitinput_folder}_{arg_postfix}', f'fitresults_asimov.hdf5')}' --result asimov --config 'utilities/styles/styles.py' --correlation --showNumbers --title CMS --subtitle WiP --titlePos 0 --scaleTextSize '0.4' -o '{os.path.join(outfolder, f'{self.fitinput_folder}_{arg_postfix}')}' --params {arg_npgamma} {arg_nplambda} 'resumTNP_gamma_nu' lumi pdf14CT18ZSymAvg pdf26CT18ZSymAvg"
    #     print(command)
    #     os.system(command)


if __name__ == "__main__":
    main()





