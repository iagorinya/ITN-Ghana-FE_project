import os
import numpy as np
from malaria.interventions.health_seeking import add_health_seeking

from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.interventions.outbreakindividual import recurring_outbreak
from malaria.reports.MalariaReport import add_summary_report

## Import custom reporters

SetupParser.default_block = 'HPC'
numseeds = 10
sim_start_year = 1960
serialize_years = 50

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=serialize_years * 365)

cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    'x_Temporary_Larval_Habitat': 0.206478,
    'Serialization_Time_Steps': [365 * serialize_years],
    'Serialization_Type': 'TIMESTEP',
    'Serialized_Population_Writing_Type': 'TIMESTEP',
    'Serialized_Population_Reading_Type': 'NONE',
    'Serialization_Mask_Node_Write': 0,
    'Serialization_Precision': 'REDUCED',
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
    'Disable_IP_Whitelist': 1,
    'x_Base_Population': 1,
    'x_Birth': 1
})
set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })


def update_cb(cb, years, serialize, ser_time_step=None):
    # Demographics
    cb.update_params({
        "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
        "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
        'Disable_IP_Whitelist': 1,
        'x_Base_Population': 1,
        'x_Birth': 1
    })

    # Report
    cb.update_params({
        'Enable_Default_Reporting': 0,
        'Enable_Demographics_Risk': 1,
        'Enable_Property_Output': 0,
        'Enable_Vector_Species_Report': 1,
        'Report_Detection_Threshold_Blood_Smear_Parasites': 50,
        "Parasite_Smear_Sensitivity": 0.02,  # 50/uL
        'RDT_Sensitivity': 0.1
    })

    return cb

# def case_management(cb, cm_cov_U5=0.30, cm_cov_adults=0.30, start_year=365):
#     # Clinical cases
#     add_health_seeking(cb, start_day=start_year,
#                        targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
#                                  'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
#                                 {'trigger': 'NewClinicalCase', 'coverage': 0.3,
#                                  'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
#                        drug=['Artemether', 'Lumefantrine'],
#                        )
#     # Severe cases
#     add_health_seeking(cb, start_day=start_year,
#                        targets=[{'trigger': 'NewSevereCase', 'coverage': 0.6,
#                                  'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
#                        drug=['Artemether', 'Lumefantrine'],
#                        broadcast_event_name='Received_Severe_Treatment')
#     return {'cm_cov_U5': cm_cov_U5,
#             'cm_cov_adults': cm_cov_adults}
#
# def itn_intervention(cb, coverage_level):
#     #itn_leak_factor = 0.9
#     seasonal_values = [0.03, 0.03, 0.01, 0.01, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.02]
#     seasonal_times = [0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 364]
#     deploy_year = 365
#     for i in [deploy_year * 54, deploy_year * 57, deploy_year * 60]:#, deploy_year * 6]:
#         add_ITN_age_season(cb,
#                            start=i,
#                            demographic_coverage=coverage_level,
#                            killing_config={
#                                "Initial_Effect": 0.7,
#                                "Decay_Time_Constant": 1460,
#                                "class": "WaningEffectExponential"},
#                            blocking_config={
#                                "Initial_Effect": 0.9,
#                                "Decay_Time_Constant": 730,
#                                "class": "WaningEffectExponential"},
#                            discard_times={
#                                "Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
#                                "Expiration_Period_Proportion_1": 0.9,
#                                "Expiration_Period_Mean_1": 365 * 1.5,
#                                "Expiration_Period_Mean_2": 3650},
#                            age_dependence={'Times': [5, 18],
#                                            'Values': [0.56, 0.2]},
#                            seasonal_dependence={"Times": seasonal_times, "Values": seasonal_values},
#                            duration=-1, birth_triggered=False)
#
#     return {'itn_coverage': coverage_level}

#recurring_outbreak(cb, start_day=180, repetitions=serialize_years)
#add_summary_report(cb, age_bins=[5, 100], start=365 * serialize_years, interval=365)

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_burnin_ITN_zone_2_{serialize_years}'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', hab_scale)]
                                for hab_scale in [0.206478]
                                for x in range(numseeds)
                                ])
# builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
#                                 ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat', hab_scale)]
#                                 for hab_scale in np.logspace(-2, np.log10(2), 7, endpoint=False)
#                                 for x in range(numseeds)
#                                 ])


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
