## Import basic python functions
import os
import pandas as pd
## Import dtk and EMOD basics functionalities
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn

## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report
from malaria.reports.MalariaReport import add_event_counter_report, add_filtered_report

## Import campaign functions
from dtk.interventions.itn import add_ITN
from malaria.interventions.health_seeking import add_health_seeking

## Required to pick up from burnin experiment
from simtools.Utilities.Experiments import retrieve_experiment

SetupParser.default_block = 'HPC'
burnin_id = "64cbc8e2-fafd-ec11-a9fa-b88303911bc1"  # UPDATE with burn-in experiment id
pull_year = 50  # year of burn-in to pick-up from
pickup_years = 10  # years of pick-up to run
numseeds = 1

SetupParser.init()
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=pickup_years * 365)

expt = retrieve_experiment(burnin_id)  # Identifies the desired burn-in experiment
# Loop through unique "tags" to distinguish between burn-in scenarios (ex. varied historical coverage levels)
ser_df = pd.DataFrame([x.tags for x in expt.simulations])
ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

cb.update_params({
    'Demographics_Filenames': [os.path.join('Namawala', 'Namawala_single_node_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Namawala', 'Namawala_single_node_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Namawala', 'Namawala_single_node_land_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Namawala', 'Namawala_single_node_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Namawala', 'Namawala_single_node_relative_humidity_daily.bin')
})

cb.update_params({
    'Serialized_Population_Reading_Type': 'READ',
    'Serialized_Population_Filenames': ['state-%05d.dtk' % (pull_year * 365)],
    # 'Serialized_Population_Path': os.path.join(ser_df["outpath"][0], 'output'),  # only use if having 1 single burnin
    'Enable_Random_Generator_From_Serialized_Population': 0,
    'Serialization_Mask_Node_Read': 0,
    'Enable_Default_Reporting': 1,
    'Disable_IP_Whitelist': 1
})

set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

"""ADDITIONAL CAMPAIGNS"""
event_list = []  ## Collect events to track in reports


# health seeeking, immediate start
def case_management(cb, cm_cov_U5=0.7, cm_cov_adults=0.5):
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': cm_cov_adults,
                                 'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'])
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': 0.85,
                                 'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'],
                       broadcast_event_name='Received_Severe_Treatment')
    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


event_list = event_list + ['Received_Treatment', 'Received_Severe_Treatment']


# ITN, start after 1 year
def itn_intervention(cb, coverage_level, day=365):
    add_ITN(cb,
            start=day,  # starts on first day of second year
            coverage_by_ages=[
                {"coverage": coverage_level, "min": 0, "max": 10},
                {"coverage": coverage_level * 0.75, "min": 10, "max": 50},
                {"coverage": coverage_level * 0.6, "min": 50, "max": 125}
            ],
            repetitions=5,  # ITN will be distributed 5 times
            tsteps_btwn_repetitions=365 * 3  # three years between ITN distributions
            )
    return {'itn_start': day,
            'itn_coverage': coverage_level}


event_list = event_list + ['Received_ITN']

"""CUSTOM REPORTS"""
add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 5, 100],
                   description='Annual_Agebin') ## (U5)

## Monthly summary report
for i in range(pickup_years):
    add_summary_report(cb, start=1 + 365 * i, interval=30,
                       duration_days=365,
                       age_bins=[0.25, 5, 120],
                       description=f'Monthly_U5_{i}')

## Enable reporters
cb.update_params({
    "Report_Event_Recorder": 0,  # No individual events tracked
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})
## Event_counter_report
add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=10000)

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management, cm_cov_U5),
                                 ModFn(itn_intervention, coverage_level=itn_cov),
                                 # Run pick-up from each unique burn-in scenario
                                 ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                                       os.path.join(row['outpath'], 'output')),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed)
                                 ]
                                for cm_cov_U5 in [0.6]
                                for itn_cov in [0, 0.9]
                                for r, row in ser_df.iterrows()
                                for seed in range(numseeds)
                                ])

user = os.getlogin()  # user initials
run_sim_args = {
    'exp_name': f'{user}_FE_2022_example_w6_test_pickup',
    'config_builder': cb,
    'exp_builder': builder
}

if __name__ == "__main__":
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
