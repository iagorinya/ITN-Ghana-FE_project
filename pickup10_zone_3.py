import os
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.interventions.itn import add_ITN
## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report, add_filtered_report
from malaria.reports.MalariaReport import add_event_counter_report
from simtools.Utilities.Experiments import retrieve_experiment
from dtk.interventions.outbreakindividual import recurring_outbreak
import pandas as pd
import numpy as np
from malaria.interventions.health_seeking import add_health_seeking
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season

SetupParser.default_block = 'HPC'
# numseeds = 3
pickup_years = 10
sim_start_year = 2010
pull_year = 50

burnin_id = 'a16c9123-f731-ed11-a9fc-b88303911bc1'
# burnin_id = '2022_08_06_21_08_47_004533'
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
    # 'x_Temporary_Larval_Habitat': 1

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
def case_management(cb, cm_cov_U5=0.30, cm_cov_adults=0.30, start_year=365):
    # Clinical cases
    for i in [start_year * 7, start_year * 8, start_year * 9]:
        add_health_seeking(cb, start_day=i,
                           targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                     'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                    {'trigger': 'NewClinicalCase', 'coverage': 0.3,
                                     'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                           drug=['Artemether', 'Lumefantrine'],
                           )
        # Severe cases
        add_health_seeking(cb, start_day=i,
                           targets=[{'trigger': 'NewSevereCase', 'coverage': 0.6,
                                     'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
                           drug=['Artemether', 'Lumefantrine'],
                           broadcast_event_name='Received_Severe_Treatment')
    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


event_list = []  ## Collect events to track in reports


def itn_intervention(cb, coverage_level=0.56):
    seasonal_values = [0.032, 0.032, 0.0378, 0.154, 0.177, 0.105, 0.25, 0.32, 0.23, 0.18, 0.032]
    seasonal_times = [0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 364]
    deploy_year = 365
    for i in [deploy_year * 6, deploy_year * 8, deploy_year * 9]:  # , deploy_year * 6, deploy_year * 9]:
        add_ITN_age_season(cb,
                           start=i,
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
                           age_dependence={'Times': [5, 18],
                                           'Values': [0.7, 0.2]},
                           seasonal_dependence={"Times": seasonal_times, "Values": seasonal_values},
                           duration=-1, birth_triggered=False)

    return {'itn_coverage': coverage_level}


event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']


# IRS, start after 1 year - single campaign
def irs_intervention(cb, coverage_level=0.60):
    deploy_year = 365
    for i in [deploy_year * 4, deploy_year * 5, deploy_year * 6, deploy_year * 7,
              deploy_year * 8, deploy_year * 9]:  # , deploy_year * 0,
        # deploy_year * 7, deploy_year * 8, deploy_year * 8, ]:
        add_IRS(cb, start=i,
                coverage_by_ages=[{"coverage": coverage_level, "min": 0, "max": 18}],
                killing_config={
                    "class": "WaningEffectBoxExponential",
                    "Box_Duration": 180,  # based on PMI data from Burkina
                    "Decay_Time_Constant": 90,  # Sumishield from Benin
                    "Initial_Effect": 0.7},
                )

    return {'irs_coverage': coverage_level}


event_list = event_list + ['Received_IRS']

# seasonal malaria chemoprevention
# def smc_intervention(cb, coverage_level=0.85, year = 365, cycles=4):
#     add_drug_campaign(cb, campaign_type='SMC', drug_code='SDX_PYR_A',
#                       coverage=coverage_level,
#                       start_days=[year*5, year *6, year*7, year * 8, year * 9],  # 4 cycles
#                       repetitions=cycles,
#                       tsteps_btwn_repetitions=30,
#                       # target_group={'agemin': 0.25, 'agemax': 5},
#                       target_group={'agemin': 0.25, 'agemax': 5},
#                       receiving_drugs_event_name='Received_SMC')
#
#     return {'smc_start': year,
#             'smc_cycle': cycles,
#             'smc_coverage': coverage_level}
#
#
# event_list = event_list + ['Received_SMC']

"""CUSTOM REPORTS"""
add_filtered_report(cb, start=0, end=pickup_years * 365)
# Summary report per agebin
add_summary_report(cb, start=0, interval=365,
                   age_bins=[0.25, 5],
                   description='U5_PfPR')

add_summary_report(cb, start=1, interval=365,
                   age_bins=[0, 5, 10, 18, 100],
                   description='Annual_Agebin')

for year in range(pickup_years):
    start_day = 0 + 365 * year
    sim_year = sim_start_year + year
    add_summary_report(cb, start=start_day, interval=30,
                       age_bins=[0.25, 5],
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
expt_name = f'{user}_FE_2022_Calibration_zone3_{serialize_year}'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management),
                                 ModFn(itn_intervention),
                                 ModFn(irs_intervention),
                                 # ModFn(smc_intervention),
                                 ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                                       os.path.join(row['outpath'], 'output')),
                                 # ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                                 #      os.path.join(ser_df[ser_df.Run_Number == seed].outpath.iloc[0], 'output')),
                                 # ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
                                 ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat',
                                       row['x_Temporary_Larval_Habitat']),
                                 ]
                                # for itn_cov in [0.3]
                                # for cm_cov_U5 in [0.2]
                                # for cm_cov_adults in [0.1]
                                # for ke in [0.06]
                                # for hab_scale in np.logspace(-2, np.log10(30), 7)
                                # for seed in range(numseeds)
                                for r, row in ser_df.iterrows()
                                ])
# (ser_df[ser_df.Run_Number==seed].
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
