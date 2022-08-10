import os
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
## Import custom reporters

SetupParser.default_block = 'HPC'
numseeds = 1
sim_start_year = 2022
serialize_years = 20

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=serialize_years * 365)
cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics_wIP.json')],
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
    'Disable_IP_Whitelist' : 1,
    'Serialization_Mask_Node_Write': 0,
    'Serialization_Precision': 'REDUCED'
})



set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_W8a'


"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x)]
                                for x in range(numseeds)])

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