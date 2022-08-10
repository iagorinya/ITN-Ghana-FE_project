## Import basic python functions
import os
## Import dtk and EMOD basics functionalities
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_larval_habitat
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
## Import custom reporters
from malaria.reports.MalariaReport import add_summary_report
from malaria.reports.MalariaReport import add_event_counter_report
## Import campaign functions
from dtk.interventions.itn import add_ITN
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.irs import add_IRS
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.interventions.malaria_vaccine import add_vaccine

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


# malaria vaccine (RTS,S), no booster start after 1 year
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


# seasonal malaria chemoprevention
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


# ITN, start after 1 year
def itn_intervention(cb, coverage_level, day=365):
    add_ITN(cb,
            start=day,  # starts on first day of second year
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
    # add_ITN_age_season(cb, start=day,
    #                    demographic_coverage=coverage_level,
    #                    killing_config={
    #                        "Initial_Effect": 0.520249973,  # LLIN Burkina
    #                        "Decay_Time_Constant": 1460,
    #                        "class": "WaningEffectExponential"},
    #                    blocking_config={
    #                        "Initial_Effect": 0.53,
    #                        "Decay_Time_Constant": 730,
    #                        "class": "WaningEffectExponential"},
    #                    discard_times={"Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
    #                                   "Expiration_Period_Proportion_1": 0.9,
    #                                   "Expiration_Period_Mean_1": 365 * 1.7,  # Burkina 1.7
    #                                   "Expiration_Period_Mean_2": 3650},
    #                    age_dependence={'Times': [0, 100],
    #                                    'Values': [0.9, 0.9]},
    #                    duration=-1, birth_triggered=False
    #                    )

    return {'itn_start': day,
            'itn_coverage': coverage_level}


event_list = event_list + ['Received_ITN']  # when using add_ITN
# event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']  # when using add_ITN_age_season


# IRS, start after 1 year - single campaign
def irs_intervention(cb, coverage_level, day=366):
    add_IRS(cb, start=day,
            coverage_by_ages=[{"coverage": coverage_level, "min": 0, "max": 100}],
            killing_config={
                "class": "WaningEffectBoxExponential",
                "Box_Duration": 180,  # based on PMI data from Burkina
                "Decay_Time_Constant": 90,  # Sumishield from Benin
                "Initial_Effect": 0.7},
            )

    return {'irs_start': day,
            'irs_coverage': coverage_level}


event_list = event_list + ['Received_IRS']

