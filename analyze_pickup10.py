from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser

from analyzer_collection import *

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials

expt_name = f'{user}_FE_2022_Calibration_zone3_50'
expt_id = 'ab1c1847-1732-ed11-a9fc-b88303911bc1'  ## change expt_id
working_dir = os.path.join('simulation_outputs')

if __name__ == "__main__":
    SetupParser.init()

    # sweep_variables = ['itn_coverage', 'Run_Number']
    sweep_variables = ['Run_Number', 'itn_coverage', 'x_Temporary_Larval_Habitat', 'cm_cov_U5',
                       'cm_cov_adults', 'irs_coverage']

    # analyzers to run
    analyzers = [
        MonthlyPfPRAnalyzerU5(expt_name=expt_name,
                              working_dir=working_dir,
                              start_year=2010,
                              end_year=2020,
                              sweep_variables=sweep_variables),
    ]
    am = AnalyzeManager(expt_id, analyzers=analyzers)
    am.analyze()
