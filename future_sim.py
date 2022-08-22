import numpy as np
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
import os
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn

## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report, add_filtered_report
from malaria.reports.MalariaReport import add_event_counter_report
from simtools.Utilities.Experiments import retrieve_experiment
from dtk.interventions.outbreakindividual import recurring_outbreak
import pandas as pd
from malaria.interventions.health_seeking import add_health_seeking

SetupParser.default_block = 'HPC'
numseeds = 5
pickup_years = 5
sim_start_year = 2020
pull_year = 50

burnin_id = '5b13fe33-7d21-ed11-a9fb-b88303911bc1'
serialize_year = 50

SetupParser.init()
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=pickup_years * 365)

expt = retrieve_experiment(burnin_id)  # Identifies the desired burn-in experiment
# Loop through unique "tags" to distinguish between burn-in scenarios (ex. varied historical coverage levels)
ser_df = pd.DataFrame([x.tags for x in expt.simulations])
ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'x_Base_Population': 1,
    'x_Birth': 1,
    'x_Temporary_Larval_Habitat': 0.14479

})

cb.update_params({
    'Serialized_Population_Reading_Type': 'READ',
    'Serialized_Population_Filenames': ['state-%05d.dtk' % (serialize_year * 365)],
    'Enable_Random_Generator_From_Serialized_Population': 0,
    'Serialization_Mask_Node_Read': 0,
    'Enable_Default_Reporting': 1
})

set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

# Add case management including two interventions
event_list = []  ## Collect events to track in reports


# health seeking, immediate start

def case_management(cb, cm_cov_U5=0.56, cm_cov_adults=0.3):
    # Clinical cases
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': 0.3,
                                 'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'],
                       )
    # Severe cases
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': cm_cov_adults,
                                 'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'],
                       broadcast_event_name='Received_Severe_Treatment')
    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


event_list = event_list + ['Received_Treatment', 'Received_Severe_Treatment']


def itn_intervention(cb, coverage_level, day=366):
    seasonal_times = [0.0, 20.0, 21.0, 30.0, 31.0, 365.0]
    seasonal_values = [0.4, 0.4, 0.5, 0.5, 0.7, 0.7]
    # for i, j in zip(seasonal_times, seasonal_values, ):
    add_ITN_age_season(cb,
                       start=day,
                       demographic_coverage=coverage_level,
                       killing_config={
                           "Initial_Effect": 0.6,
                           "Decay_Time_Constant": 1460,
                           "class": "WaningEffectExponential"},
                       blocking_config={
                           "Initial_Effect": 0.9,
                           "Decay_Time_Constant": 730,
                           "class": "WaningEffectExponential"},
                       discard_times={
                           "Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
                           "Expiration_Period_Proportion_1": 0.9,
                           "Expiration_Period_Mean_1": 365 * 1.5,
                           "Expiration_Period_Mean_2": 3650},
                       age_dependence={'Times': [0, 5, 18],
                                       'Values': [1, 0.7, 0.2]},
                       seasonal_dependence={"Times": seasonal_times, "Values": seasonal_values},
                       duration=-1, birth_triggered=False),

    return {'itn_start': day,
            'itn_coverage': coverage_level}


event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']


# IRS, start after 1 year - single campaign
def irs_intervention(cb, coverage_level, day=366):
    add_IRS(cb, start=day,
            coverage_by_ages=[{"coverage": coverage_level, "min": 0, "max": 100}],
            killing_config={
                "class": "WaningEffectBoxExponential",
                "Box_Duration": 180,  # based on PMI data from Burkina
                "Decay_Time_Constant": 90,  # Sumishield from Benin
                "Initial_Effect": 0.7},
            )

    return {'irs_start': day,
            'irs_coverage': coverage_level}


event_list = event_list + ['Received_IRS']
# # IRS intervention is only added for the Savannah Ecological zone starting 2014 (Onces a year at one year interval)
# def irs_intervention(cb, KE=0.15):
#     day = [1464, 1830, 2196, 2562]
#     coverage_level = [0.06, 0.06, 0.06, 0.06]
#     KEs = [KE] * 5
#     for i, j, k in zip(day, coverage_level, KEs):
#         add_IRS(cb, start=i,  # starts on first day of second year
#                 coverage_by_ages=[
#                     {"coverage": coverage_level, "min": 0, "max": 100}],
#                 killing_config={
#                     "class": "WaningEffectBoxExponential",
#                     "Box_Duration": 120,
#                     "Decay_Time_Constant": 120,
#                     "Initial_Effect": k})
#     return {'Killing_Effect': KE,
#             'irs_coverage': coverage_level}
#
#
# event_list = event_list + ['Received_IRS']

"""CUSTOM REPORTS"""
add_filtered_report(cb, start=0, end=pickup_years * 365)
# Summary report per agebin
add_summary_report(cb, start=0, interval=365,
                   age_bins=[0.25, 5],
                   description='U5_PfPR')

add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 5, 10, 18],
                   description='Annual_Agebin')

for year in range(pickup_years):
    start_day = 0 + 365 * year
    sim_year = sim_start_year + year
    add_summary_report(cb, start=start_day, interval=30,
                       age_bins=[0, 125],
                       description=f'Monthly_U5_{sim_year}')

# Enable reporters
cb.update_params({
    "Report_Event_Recorder": 1,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})
# Event_counter_report
add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=10000)
recurring_outbreak(cb, start_day=180, repetitions=pickup_years)

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_pickup_ITN_calibration_9_{serialize_year}'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management),# cm_cov_U5, cm_cov_adults),
                                 ModFn(itn_intervention, coverage_level=itn_cov),
                                 ModFn(irs_intervention, coverage_level=irs_cov),
                                 ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                                       os.path.join(row['outpath'], 'output')),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
                                 ModFn(DTKConfigBuilder.set_param, 'Scenario', 'Basic')  # optional

                                 ]
                                #for cm_cov_U5 in [0.56]
                                #for cm_cov_adults in [0.2]
                                for itn_cov in [0, 0.4, 0.56, 0.8]
                                for irs_cov in [0, 0.4, 0.56, 0.8]
                                for seed in range(numseeds)
                                for r, row in ser_df.iterrows()
                                ])

run_sim_args = {
    'exp_name': expt_name,
    'config_builder': cb,
    'exp_builder': builder
}
# If you prefer running with `python example_sim.py`, you will need the following block
if __name__ == "__main__":
    # SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
