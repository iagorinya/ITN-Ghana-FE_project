import os
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.interventions.itn import add_ITN
## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report
from malaria.reports.MalariaReport import add_event_counter_report
from simtools.Utilities.Experiments import retrieve_experiment
import pandas as pd
import numpy as np
from malaria.interventions.health_seeking import add_health_seeking

SetupParser.default_block = 'HPC'
numseeds = 2
pickup_years = 11
sim_start_year = 2010
pull_year = 50

burnin_id = '8ea49fc8-9b1e-ed11-a9fb-b88303911bc1'
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
    'x_Birth': 1
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
event_list = []  ## Collect events to track in reports


# health seeking, immediate start

def case_management(cb, cm_cov_U5=0.7, cm_cov_adults=0.5):
    # Clinical cases
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': 0.3,
                                 'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'])
    # Severe cases
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': 0.49,
                                 'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'],
                       broadcast_event_name='Received_Severe_Treatment')
    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


event_list = event_list + ['Received_Treatment', 'Received_Severe_Treatment']


# Add intervention
def itn_intervention(cb, coverage_level):
    add_ITN(cb,
            start=366,  # starts on first day of second year
            coverage_by_ages=[
                {"coverage": coverage_level, "min": 0, "max": 5},  # Highest coverage for 0-10 years old
                {"coverage": coverage_level * 0.75, "min": 5, "max": 100},
                # 25% lower than for children for 10-50 years old
                # {"coverage": coverage_level * 0.6, "min": 50, "max": 125}
                # 40% lower than for children for everyone else
            ],
            repetitions=5,  # ITN will be distributed 5 times
            tsteps_btwn_repetitions=365 * 3  # three years between ITN distributions
            )
    return {'itn_coverage': coverage_level}


event_list = event_list + ['Received_ITN']

"""CUSTOM REPORTS"""
# add_filtered_report(cb, start=0, end=years * 365)
# Summary report per agebin
add_summary_report(cb, start=366, interval=365,
                   age_bins=[5, 100],
                   description='U5_PfPR')

add_summary_report(cb, start=366, interval=365,
                   age_bins=[0.25, 5, 100],
                   description='Annual_U5')

# add_summary_report(cb, start=1, interval=365,
#                    age_bins=[0.25, 10, 100],
#                    description='Annual_U10')

for year in range(pickup_years):
    start_day = 0 + 365 * year
    sim_year = sim_start_year + year
    add_summary_report(cb, start=start_day, interval=30,
                       age_bins=[0.25, 2, 5, 10, 15, 20, 100, 120],
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

# run_sim_args is what the `dtk run` command will look for
user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_pickup_ITN_calibration_4_{serialize_year}'

"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(case_management, cm_cov_U5),
                                 ModFn(itn_intervention, coverage_level=itn_cov),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                                       os.path.join(row['outpath'], 'output')),
                                 ModFn(DTKConfigBuilder.set_param, 'x_Temporary_Larval_Habitat',
                                       row['x_Temporary_Larval_Habitat']),
                                 ]
                                for itn_cov in [0.69]
                                for cm_cov_U5 in [0.5]
                                # for hab_scale in np.logspace(-2, np.log10(30), 50)
                                for x in range(numseeds)
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
