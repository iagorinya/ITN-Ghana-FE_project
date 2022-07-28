import os

from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report
from malaria.reports.MalariaReport import add_event_counter_report
##Import campaign functions
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from dtk.interventions.itn import add_ITN
from Solution_scripts.run_exampleSim_w3b import itn_intervention, irs_intervention

SetupParser.default_block = 'HPC'
numseeds = 3
years = 3
event_list = []  ## Collect events to track in reports
sim_start_year = 2022

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=years * 365)
cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'x_Base_Population': 1,
    'x_Birth': 1,
    'x_Temporary_Larval_Habitat': 1,

})

set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

# Add case management including two interventions
event_list = []  ## Collect events to track in reports

# Clinical cases
def case_management(cb, cm_cov_U5=0.7, cm_cov_adults=0.5):
    add_health_seeking(cb, start_day=366,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': cm_cov_adults,
                                 'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'])
    # Severe cases
    add_health_seeking(cb, start_day=366,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': 0.49,
                                 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'],
                       broadcast_event_name='Received_Severe_Treatment')
    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}

event_list = event_list + ['Received_Treatment', 'Received_Severe_Treatment']

# Adding ITNs
# ITN, start after 1 year
"""add_ITN_age_season"""

def ITN_management(cb, start_day=185, coverage_level=0.7):
    add_ITN_age_season(cb, start=start_day,
                       demographic_coverage=coverage_level,
                       killing_config={
                           "Initial_Effect": 0.520249973,  # LLIN Burkina
                           "Decay_Time_Constant": 1460,
                           "class": "WaningEffectExponential"},
                       blocking_config={
                           "Initial_Effect": 0.53,
                           "Decay_Time_Constant": 730,
                           "class": "WaningEffectExponential"},
                       discard_times={"Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
                                      "Expiration_Period_Proportion_1": 0.9,
                                      "Expiration_Period_Mean_1": 365 * 1.7,  # Burkina 1.7
                                      "Expiration_Period_Mean_2": 3650},
                       age_dependence={'Times': [0, 100],
                                       'Values': [0.9, 0.9]},
                       duration=-1, birth_triggered=False
                       )
event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']  # when using add_ITN_age_season

"""CUSTOM REPORTS"""
# add_filtered_report(cb, start=0, end=years * 365)
## Summary report per agebin
add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 2, 5, 10, 15, 20, 100, 120],
                   description='Annual_Agebin')

add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 5, 100],
                   description='Annual_U5')

add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 10, 100],
                   description='Annual_U10')

for year in range(years):
    start_day = 365 + 365 * year
    sim_year = sim_start_year + year
     add_summary_report(cb, start=start_day, interval=30,
                       age_bins=[0.25, 5, 100],
                       description=f'Monthly_U5_{sim_year}')

   ## Enable reporters
cb.update_params({
    "Report_Event_Recorder": 1,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})
# Event_counter_report
add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=10000)

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_ITN_Ghana'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management, cm_cov_U5),
                                 ModFn(itn_intervention, coverage_level=itn_cov),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', x)
                                 ]
                                for cm_cov_U5 in [0.5, 0.6 ]
                                for itn_cov in [0.4]
                                for x in range(numseeds)
                                ])

run_sim_args = {
    'exp_name': expt_name,
    'config_builder': cb,
    'exp_builder': builder
}

# If you prefer running with `python example_sim.py`, you will need the following block
if __name__ == "__main__":
    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
