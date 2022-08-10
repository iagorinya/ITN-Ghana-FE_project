import os
import numpy as np
import pandas as pd
from dtk.interventions.itn import add_ITN
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from malaria.reports.MalariaReport import add_event_counter_report, add_summary_report
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn

# This block will be used unless overridden on the command-line
from simtools.Utilities.Experiments import retrieve_experiment

SetupParser.default_block = 'HPC'
burnin_id = "059e19a7-3e17-ed11-a9fb-b88303911bc1"  # UPDATE with burn-in experiment id
pull_year = 5  # year of burn-in to pick-up from
pickup_years = 10  # years of pick-up to run
numseeds = 1

SetupParser.init()
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=pickup_years * 365)

expt = retrieve_experiment(burnin_id)  # Identifies the desired burn-in experiment
# Loop through unique "tags" to distinguish between burn-in scenarios (ex. varied historical coverage levels)
ser_df = pd.DataFrame([x.tags for x in expt.simulations])
ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

"""Demographics"""
cb.update_params({
    'Demographics_Filenames': [os.path.join('Mozambique', 'Mozambique_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Mozambique', 'Malema_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Mozambique', 'Malema_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Mozambique', 'Malema_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Mozambique', 'Malema_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE'
})

"""Serialization"""
cb.update_params({
    'Serialized_Population_Reading_Type': 'READ',
    'Serialized_Population_Filenames': ['state-%05d.dtk' % (pull_year * 365)],
    'Enable_Random_Generator_From_Serialized_Population': 0,
    'Serialization_Mask_Node_Read': 0,
    'Enable_Default_Reporting': 1
})

set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })


# Add intervention
def itn_intervention(cb, coverage_level):
    add_ITN(cb,
            start=1,  # starts on first day of second year
            coverage_by_ages=[
                {"coverage": coverage_level, "min": 0, "max": 10},  # Highest coverage for 0-10 years old
                {"coverage": coverage_level * 0.75, "min": 10, "max": 50},
                # 25% lower than for children for 10-50 years old
                {"coverage": coverage_level * 0.6, "min": 50, "max": 125}
                # 40% lower than for children for everyone else
            ],
            repetitions=5,  # ITN will be distributed 5 times
            tsteps_btwn_repetitions=365 * 3  # three years between ITN distributions
            )
    return {'itn_coverage': coverage_level}


event_list = []
event_list = event_list + ['Received_ITN']
cb.update_params({
    "Report_Event_Recorder": 1,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})

# Report_Event_Counter
# add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=pickup_years*365)

for i in range(pickup_years):
    add_summary_report(cb, start=1+365*i, interval=30,
                       duration_days=365,
                       age_bins=[5, 100],
                       # age_bins=[0.25, 5, 120],
                       description=f'Monthly_U5_{2010+i}')

# # run_sim_args is what the `dtk run` command will look for
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', m)]
                                for m in np.logspace(1, np.log10(30), 5)
                                for n in range(numseeds)])

user = os.getlogin()  # user initials
run_sim_args = {
    'exp_name': f'{user}_FE_2022_pickupMoza4',
    'config_builder': cb,
    'exp_builder': builder
}

# If you prefer running with `python example_sim.py`, you will need the following block
if __name__ == "__main__":
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
