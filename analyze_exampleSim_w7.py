from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser

from analyzer_collection import *

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials
# expt_name = f'{user}_FE_2022_example_w7b'
# expt_id = 'c9eccdcd-b711-ed11-a9fb-b88303911bc1'  ## change expt_id
# working_dir = os.path.join('simulation_outputs')
expt_name = f'{user}_FE_2022_w7pick_up50'  ## change w6a or w6b
expt_id = '22298e88-2018-ed11-a9fb-b88303911bc1'  ## change expt_id
working_dir = os.path.join('simulation_outputs')

serialize_years = 50  # Same as in run_exampleSim_w6a.py
step = 'pickup'  # 'pi


if __name__ == "__main__":
    SetupParser.init()

    sweep_variables = ['itn_coverage', 'Run_Number']

    # analyzers to run
    analyzers = [
        MonthlyPfPRAnalyzerU5(expt_name=expt_name,
                              working_dir=working_dir,
                              start_year=2020,
                              end_year=2030,
                              sweep_variables=sweep_variables),
    ]
    am = AnalyzeManager(expt_id, analyzers=analyzers)
    am.analyze()
