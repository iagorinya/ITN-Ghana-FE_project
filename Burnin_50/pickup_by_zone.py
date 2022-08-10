import os
from dtk.interventions.outbreakindividual import recurring_outbreak
import pandas as pd
import numpy as np
from dtk.interventions.itn import add_ITN
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from malaria.reports.MalariaReport import add_event_counter_report, add_summary_report
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.Utilities.Experiments import retrieve_experiment
## Import custom reporters

SetupParser.default_block = 'HPC'
sim_start_year = 2010  # for forward simulations with future intervention scenarios
numseeds = 1
pickup_years = 10  # 2010 to 2019
burnin_id = "7bcae914-2f15-ed11-a9fb-b88303911bc1"  # UPDATE with burn-in experiment id

SetupParser.init()
expt = retrieve_experiment(burnin_id)  # Identifies the desired burn-in experiment
# Loop through unique "tags" to distinguish between burn-in scenarios (ex. varied historical coverage levels)
ser_df = pd.DataFrame([x.tags for x in expt.simulations])
ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=pickup_years * 365)
cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'x_Base_Population': 1,
    'x_Birth': 1,
    'x_Temporary_Larval_Habitat': 1
})
set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

def update_cb(cb, years, serialize, ser_time_step = None):

    # # Logging
    # cb.update_params({
    #     'logLevel_JsonConfigurable' : 'ERROR',
    #     'logLevel_VectorHabitat' : 'ERROR',
    #     'logLevel_StandardEventCoordinator' : 'ERROR',
    #     'logLevel_SusceptibilityMalaria' : 'ERROR'
    # })

    # Demographics
    cb.update_params({
        "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
        "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
        'Disable_IP_Whitelist' : 1,
        'x_Base_Population': 1,
        'x_Birth': 1
    })

    # Serialization
    if ser_time_step is None:
        ser_time_step = [365*years]
    if serialize :
        cb.update_params({
            'Serialization_Time_Steps' : ser_time_step,
            'Serialization_Type': 'TIMESTEP',
            'Serialized_Population_Writing_Type': 'TIMESTEP',
            'Serialization_Mask_Node_Write': 0,
            # 0 corresponds to the previous version default: the same larval habitat
            # parameters will be used in the burnin and pickup (from the burnin config)
            'Serialization_Precision': 'REDUCED'
        })
    else:
        cb.update_params({
            'Serialization_Type': 'NONE',
            'Serialized_Population_Writing_Type': 'NONE'
        })

    # Report
    cb.update_params( {
        'Enable_Default_Reporting': 0,
        'Enable_Demographics_Risk': 1,
        'Enable_Property_Output': 0,
        'Enable_Vector_Species_Report': 1,
        'Report_Detection_Threshold_Blood_Smear_Parasites': 50,
        "Parasite_Smear_Sensitivity": 0.02,  # 50/uL
        'RDT_Sensitivity': 0.1
    })

    return cb
#
# def set_vectors_and_habitats_sweep(cb, larv_hab_multiplier):
#     # Vector
#     cb.update_params({
#         "Vector_Species_Names": ['arabiensis', 'funestus', 'gambiae'],
#         'x_Temporary_Larval_Habitat': larv_hab_multiplier
#     })
#     set_species(cb, ["arabiensis", "funestus", "gambiae"])
#     set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
#                             "funestus": {'WATER_VEGETATION': 4e8},
#                             "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
#                             })
#     return {'larv_hab_multiplier': larv_hab_multiplier}
#
# def set_input_files(cb, my_ds, demographic_suffix='_2.5arcmin', climate_suffix='_30arcsec_air'):
#
#     if demographic_suffix is not None:
#         if not demographic_suffix.startswith('_') and not demographic_suffix == '':
#             demographic_suffix = '_' + demographic_suffix
#
#     if climate_suffix is not None:
#         if not climate_suffix.startswith('_') and not climate_suffix == '':
#             climate_suffix = '_' + climate_suffix
#
#     if demographic_suffix is not None:
#         cb.update_params({
#             'Demographics_Filenames': [os.path.join(my_ds, f'{my_ds}{demographic_suffix}_demographics%s.json')]
#         })
#     if climate_suffix is not None:
#         cb.update_params({
#             "Air_Temperature_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_air_temperature_daily%s.bin' ),
#             "Land_Temperature_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_air_temperature_daily%s.bin' ),
#             "Rainfall_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_rainfall_daily%s.bin' ),
#             "Relative_Humidity_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_relative_humidity_daily%s.bin')
#         })
#
#     return {'DS_Name': my_ds}

recurring_outbreak(cb, start_day=180, repetitions=pickup_years)
add_summary_report(cb, age_bins=[5, 100], start=365*(pickup_years-1), interval=365)

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
# add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=pull_year*365)

for i in range(pickup_years):
    add_summary_report(cb, start=1+365*i, interval=30,
                       duration_days=365,
                       age_bins=[0.25, 5, 120],
                       description=f'Habitat_multiplier{2010+i}')

# run_sim_args is what the `dtk run` command will look for
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', hab_scale)]
                                for hab_scale in np.logspace(-2, np.log10(30), 7)
                                for x in range(numseeds)
                                ])

user = os.getlogin()  # user initials
run_sim_args = {
    'exp_name': f'{user}_FE_2022_Pickup{pickup_years}',
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