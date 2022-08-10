import os
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
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.irs import add_IRS
from dtk.interventions.novel_vector_control import add_larvicides
from malaria.interventions.malaria_vaccine import add_vaccine

from Solution_scripts.run_exampleSim_w3b import itn_intervention, irs_intervention

SetupParser.default_block = 'LOCAL'
numseeds = 3
years = 3
sim_start_year = 2022

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=years * 365)
cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'x_Base_Population': 2,
    'x_Birth': 1,
    'x_Temporary_Larval_Habitat': 1
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
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': 0.7,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': 0.5,
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

# malaria vaccine (RTS,S), as mass campaign at 80% at a single day (hypothetical example only!)
def rtss_intervention(cb, coverage_level, day=366, agemin=274, agemax=275, initial_efficacy=0.8):
    add_vaccine(cb,
                vaccine_type='RTSS',
                vaccine_params={"Waning_Config":
                                    {"Initial_Effect": initial_efficacy,
                                     "Decay_Time_Constant": 592.4066512,
                                     "class": 'WaningEffectExponential'}},
                start_days=[day],
                coverage=coverage_level,
                repetitions=1,
                tsteps_btwn_repetitions=-1,
                target_group={'agemin': agemin, 'agemax': agemax})  # children 9 months of age

    return {'rtss_start': day,
            'rtss_coverage': coverage_level,
            'rtss_initial_effect': initial_efficacy}

event_list = event_list + ['Received_Vaccine']

# seasonal malaria chemoprevention, start after 1 year
def smc_intervention(cb, coverage_level, day=366, cycles=4):
    add_drug_campaign(cb, campaign_type='SMC', drug_code='SPA',
                      coverage=coverage_level,
                      start_days=[day],
                      repetitions=cycles,
                      tsteps_btwn_repetitions=30,
                      target_group={'agemin': 0.25, 'agemax': 5},
                      receiving_drugs_event_name='Received_SMC')

    return {'smc_start': day,
            'smc_coverage': coverage_level}

event_list = event_list + ['Received_SMC']

# Adding ITNs
# ITN, start after 1 year
"""Select either add_ITN or add_ITN_age_season"""
add_ITN(cb,
        start=366,  # starts on first day of second year
        coverage_by_ages=[
            {"coverage": 1, "min": 0, "max": 10},  # 100% for 0-10 years old
            {"coverage": 0.75, "min": 10, "max": 50},  # 75% for 10-50 years old
            {"coverage": 0.6, "min": 50, "max": 125}  # 60% for everyone else
        ],
        repetitions=5,  # ITN will be distributed 5 times
        tsteps_btwn_repetitions=366 * 3  # assume annual ITN distributions instead of 1 year in example
        )
event_list = event_list + ['Received_ITN']
#
# IRS, start after 1 year - single campaign
add_IRS(cb,
        start=366,  # IRS occurs on first day of second year
        coverage_by_ages=[
            {"coverage": 1, "min": 0, "max": 10},  # 100% for 0-10 years old
            {"coverage": 0.75, "min": 11, "max": 50},  # 75% for 11-50 years old
            {"coverage": 0.6, "min": 51, "max": 125}  # 60% for everyone else
        ],
        killing_config={
            "class": "WaningEffectBoxExponential",
            "Box_Duration": 60,
            "Decay_Time_Constant": 120,
            "Initial_Effect": 0.6
        }
        )
event_list = event_list + ['Received_IRS']
# #
# Larviciding, start after 1 year - single campaign
add_larvicides(cb, start_day=366,
               habitat_target='CONSTANT',
               killing_initial=0.4,
               killing_decay=150)

"""CUSTOM REPORTS"""
# add_filtered_report(cb, start=0, end=years * 365)
## Summary report per agebin
add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 2, 5, 10, 15, 20, 100, 120],
                   description='Annual_Agebin')

## Enable reporters
cb.update_params({
    "Report_Event_Recorder": 0,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})
# Event_counter_report
add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=10000)

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_example_w3b'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management, cm_cov_U5),
                                 ModFn(smc_intervention, coverage_level=smc_cov),
                                 ModFn(rtss_intervention, coverage_level=rtss_cov),
                                 ModFn(itn_intervention, coverage_level=itn_cov),
                                 ModFn(irs_intervention, coverage_level=irs_cov),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', x)
                                 ]
                                for cm_cov_U5 in [0.6]
                                for smc_cov in [0, 1]
                                for rtss_cov in [0]
                                for itn_cov in [0]
                                for irs_cov in [0]
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
