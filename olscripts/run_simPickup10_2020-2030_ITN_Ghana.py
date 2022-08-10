## Import basic python functions
import os
import pandas as pd
import numpy as np
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
from dtk.interventions.itn_age_season import add_ITN_age_season
from malaria.interventions.health_seeking import add_health_seeking

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'
sim_start_year = 2020  # for forward simulations with future intervention scenarios
numseeds = 3
pickup_years = 2  # 2020 to 2030

user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_pick_up{pickup_years}'  # give simulation a systematic, meaningful name

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM', Simulation_Duration=pickup_years * 365)

cb.update_params({
    'Demographics_Filenames': [os.path.join('Ghana', 'Ghana_2.5arcmin_demographics.json')],
    "Air_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Land_Temperature_Filename": os.path.join('Ghana', 'Ghana_30arcsec_air_temperature_daily.bin'),
    "Rainfall_Filename": os.path.join('Ghana', 'Ghana_30arcsec_rainfall_daily.bin'),
    "Relative_Humidity_Filename": os.path.join('Ghana', 'Ghana_30arcsec_relative_humidity_daily.bin'),
    "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
    'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE'
})

set_species(cb, ["arabiensis", "funestus", "gambiae"])
set_larval_habitat(cb, {"arabiensis": {'TEMPORARY_RAINFALL': 7.5e9, 'CONSTANT': 1e7},
                        "funestus": {'WATER_VEGETATION': 4e8},
                        "gambiae": {'TEMPORARY_RAINFALL': 8.3e8, 'CONSTANT': 1e7}
                        })

"""ADDITIONAL CAMPAIGNS"""

# health seeeking, immediate start
def add_ds_hs(cb, cm_df, my_ds):
    df = cm_df[cm_df['repDS'] == my_ds]

    # Multiple start days and coverage levels per setting if case management were to increase,
    # here assumings single values for constant case management over time
    cm_cov_U5 = df['U5_coverage'][0]
    cm_cov_adults = df['adult_coverage'][0]
    cm_start = df['simday'][0]

    add_health_seeking(cb, start_day=cm_start,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5,
                                 'agemin': 0, 'agemax': 5, 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': cm_cov_adults,
                                 'agemin': 5, 'agemax': 100, 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'])

    add_health_seeking(cb, start_day=cm_start,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': 0.85,
                                 'agemin': 0, 'agemax': 100, 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'],
                       broadcast_event_name='Received_Severe_Treatment')
    return {'repDS': my_ds,
            'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


def adjust_itn_seasonals(simday):
    ##TODO hardcoded, but would also be setting specific
    itn_seasonal_values = [0.032, 0.032, 0.0378, 0.154, 0.177, 0.105, 0.25, 0.32, 0.23, 0.18, 0.032]
    itn_seasonal_months = [0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 364]

    seasonal_scales = [x / max(itn_seasonal_values) for x in itn_seasonal_values]
    seasonal_offset = simday % 365
    seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in itn_seasonal_months]

    zipped_lists = zip(seasonal_times, seasonal_scales)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
    if seasonal_times[0] > 0:
        seasonal_times.insert(0, 0)
        seasonal_scales.insert(0, seasonal_scales[-1])

    return (seasonal_times, seasonal_scales)


def ITN_age_season_template(cb, start, demographic_coverage,
                            killing_rate, blocking_rate,
                            age_bin, usages, duration=-1,
                            itn_retention_in_yr=1.51,
                            birth_triggered=False,
                            trigger_condition_list=None,
                            ind_property_restrictions=None):
    seasonal_times, seasonal_scales = adjust_itn_seasonals(start)

    add_ITN_age_season(cb, start=start,
                       demographic_coverage=demographic_coverage,
                       killing_config={
                           "Initial_Effect": killing_rate,
                           "Decay_Time_Constant": 1460,
                           "class": "WaningEffectExponential"},
                       blocking_config={
                           "Initial_Effect": blocking_rate,
                           "Decay_Time_Constant": 730,
                           "class": "WaningEffectExponential"},
                       discard_times={"Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
                                      "Expiration_Period_Proportion_1": 0.9,
                                      "Expiration_Period_Mean_1": 365 * itn_retention_in_yr,
                                      "Expiration_Period_Mean_2": 3650},
                       age_dependence={'Times': age_bin,
                                       'Values': usages},
                       seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales},
                       duration=duration, birth_triggered=birth_triggered,
                       trigger_condition_list=trigger_condition_list,
                       ind_property_restrictions=ind_property_restrictions
                       )


def add_itn_by_row(cb, row, itn_cov_sweep=None):
    itn_cov_cols = ['U5_ITN_use', 'six_nine_ITN_use', 'ten_eighteen_ITN_use', 'over_eighteen_ITN_use']
    itn_cov_age_bin = [0, 5, 10, 18]
    itn_leak_factor = 0.9

    usages = [row[x] for x in itn_cov_cols]
    ## itn_cov_sweep  was added to allow optional sweeping through itn coverage, assuming the same coverages for all years!
    ## Otherwise a scaling factor, could be used and applied on the rows of itn_df
    if itn_cov_sweep is not None:
        coverage_all = itn_cov_sweep
    else:
        coverage_all = np.max(usages)
    if coverage_all == 0:
        coverage_all = 1
    usages = [x / coverage_all for x in usages]

    ITN_age_season_template(cb, start=row['simday'],
                            demographic_coverage=coverage_all,
                            killing_rate=row['kill_rate'],
                            blocking_rate=row['blocking_rate'],
                            age_bin=itn_cov_age_bin,
                            usages=[x * itn_leak_factor for x in usages])
    return {'itn_coverage_max': coverage_all}


def add_itn_antenatal_by_row(cb, row):
    itn_leak_factor = 0.9
    # for r, row in itn_anc_df.iterrows() :
    ITN_age_season_template(cb, start=row['simday'],
                            demographic_coverage=row['coverage'],
                            killing_rate=row['kill_rate'],
                            blocking_rate=row['blocking_rate'],
                            age_bin=[0, 5],
                            usages=[itn_leak_factor for i in range(2)],
                            birth_triggered=True)


"""REGULAR ITN MASS CAMPAIGNS"""


def add_ds_itns(cb, itn_df, my_ds, itn_cov_sweep=None):
    df = itn_df[itn_df['repDS'].str.upper() == my_ds.upper()]
    df = df.drop_duplicates()
    for r, row in df.iterrows():
        if 'blocking_rate' not in row.index:
            row['blocking_rate'] = row['block_initial']
        add_itn_by_row(cb, row, itn_cov_sweep=itn_cov_sweep)

    return {'itn_cov_sweep': itn_cov_sweep}  ## dummy function requires dictionary output


"""ADDITIONAL NETS"""

def add_ds_itns_addtnl(cb, itn_addtnl_df, my_ds):
    df = itn_addtnl_df[itn_addtnl_df['repDS'].str.upper() == my_ds.upper()]
    df = df.drop_duplicates()
    nets = len(df)
    for r, row in df.iterrows():
        itn_type = row['type']
        if 'blocking_rate' not in row.index:
            row['blocking_rate'] = row['block_initial']
        #if itn_type == 'pregnant':
            #add_itn_pregnant_by_row(cb, row)
        if itn_type == 'antenatal':
            add_itn_antenatal_by_row(cb, row)
        else:
            raise ValueError('Specify itn_type, allowed values are pregnant or antenatal')

    return {'itn_routine': 'itn_routine'}  ## dummy function requires dictionary output



"""CUSTOM REPORTS"""
add_filtered_report(cb, start=0, end=pickup_years * 365)
## Summary report for 3month to 5 years (U5)
add_summary_report(cb, start=1, interval=365,
                   age_bins=[0.25, 5, 100],
                   description='Annual_U5')

for year in range(pickup_years):
    start_day = 365 + 365 * year
    sim_year = sim_start_year + year
    add_summary_report(cb, start=start_day, interval=30,
                       age_bins=[0.25, 5, 100],
                       description=f'Monthly_U5_{sim_year}')

## Enable reporters
event_list = ['Received_Treatment', 'Received_Severe_Treatment']
event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']

cb.update_params({
    "Report_Event_Recorder": 1,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Ignore_Events_In_List": 0,
    "Report_Event_Recorder_Events": event_list,
    'Custom_Individual_Events': event_list
})

## Event_counter_report
add_event_counter_report(cb, event_trigger_list=event_list, start=0, duration=10000)

"""Load example  dataframes for setting specific input values"""
repDS_list = ['GHA_ecozone1', 'GHA_ecozone2', 'GHA_ecozone3']
inputs_path = os.path.join('./', 'input/Ghana')  # TODO add custom path to input files
cm_df = pd.read_csv(os.path.join(inputs_path, 'cm_scenario1_example_2020_2030.csv'))
itn_df = pd.read_csv(os.path.join(inputs_path, 'itn_scenario1_example_2020_2030.csv'))
itn_anc_df = pd.read_csv(os.path.join(inputs_path,'itn_anc_scenario1_example_2020_2030.csv'))  ##


"""BUILDER"""
builder = ModBuilder.from_list([[ModFn(add_ds_hs, cm_df, my_ds),
                                 ModFn(add_ds_itns, itn_df, my_ds, itn_cov_sweep),  # if itn_cov_sweep is not provided, default setting specific values are used
                                 ModFn(add_ds_itns_addtnl, itn_anc_df, my_ds),
                                 # Run pick-up from each unique burn-in scenario
                                 # ModFn(DTKConfigBuilder.set_param,
                                 #      'Serialized_Population_Path',
                                 #      os.path.join(ser_df[ser_df.repDS == repDS].outpath.iloc[0], 'output')),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed)
                                 ]
                                for itn_cov_sweep in [0, 0.8]
                                # optional, allows custom sweep and variation in total ITN coverage
                                for my_ds in [repDS_list[0]]  # for testing run for 1 setting only
                                for seed in range(numseeds)
                                ])

# run_sim_args is what the `dtk run` command will look for
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
