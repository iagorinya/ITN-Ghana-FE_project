import os
import numpy as np
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.interventions.outbreakindividual import recurring_outbreak
from malaria.reports.MalariaReport import add_summary_report
## Import custom reporters

SetupParser.default_block = 'HPC'
numseeds = 1
sim_start_year = 1960
serialize_years = 50

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=serialize_years * 365)
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
    'Serialization_Time_Steps': [365 * serialize_years],
    'Serialization_Type': 'TIMESTEP',
    'Serialized_Population_Writing_Type': 'TIMESTEP',
    'Serialized_Population_Reading_Type': 'NONE',
    'Serialization_Mask_Node_Write': 0,
    'Serialization_Precision': 'REDUCED'
})
set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

def update_cb(cb, years, serialize, ser_time_step = None):

    # Logging
    cb.update_params({
        'logLevel_JsonConfigurable' : 'ERROR',
        'logLevel_VectorHabitat' : 'ERROR',
        'logLevel_StandardEventCoordinator' : 'ERROR',
        'logLevel_SusceptibilityMalaria' : 'ERROR'
    })

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

def set_vectors_and_habitats_sweep(cb, larv_hab_multiplier):
    # Vector
    cb.update_params({
        "Vector_Species_Names": ['arabiensis', 'funestus', 'gambiae'],
        'x_Temporary_Larval_Habitat': larv_hab_multiplier
    })
    set_species(cb, ["arabiensis", "funestus", "gambiae"])
    set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                            "funestus": {'WATER_VEGETATION': 4e8},
                            "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                            })
    return {'larv_hab_multiplier': larv_hab_multiplier}

def set_input_files(cb, my_ds, demographic_suffix='_2.5arcmin', climate_suffix='_30arcsec_air'):

    if demographic_suffix is not None:
        if not demographic_suffix.startswith('_') and not demographic_suffix == '':
            demographic_suffix = '_' + demographic_suffix

    if climate_suffix is not None:
        if not climate_suffix.startswith('_') and not climate_suffix == '':
            climate_suffix = '_' + climate_suffix

    if demographic_suffix is not None:
        cb.update_params({
            'Demographics_Filenames': [os.path.join(my_ds, f'{my_ds}{demographic_suffix}_demographics%s.json')]
        })
    if climate_suffix is not None:
        cb.update_params({
            "Air_Temperature_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_air_temperature_daily%s.bin' ),
            "Land_Temperature_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_air_temperature_daily%s.bin' ),
            "Rainfall_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_rainfall_daily%s.bin' ),
            "Relative_Humidity_Filename": os.path.join(my_ds, f'{my_ds}{climate_suffix}_relative_humidity_daily%s.bin')
        })

    return {'DS_Name': my_ds}

recurring_outbreak(cb, start_day=180, repetitions=serialize_years)
add_summary_report(cb, age_bins=[5, 100], start=365*(serialize_years-1), interval=365)

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_burnin{serialize_years}'


"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', hab_scale)]
                                for hab_scale in np.logspace(-2, np.log10(30), 7)
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