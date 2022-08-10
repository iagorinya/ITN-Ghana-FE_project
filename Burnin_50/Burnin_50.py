from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
import numpy as np
import os

from analyzer_collection import *

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials
expt_name = f'{user}_test3'
expt_id = '0ac5096d-ea10-ed11-a9fb-b88303911bc1'  ## change expt_id
working_dir = os.path.join('simulation_outputs')


def plot_inset_chart(channels_inset_chart, sweep_variables) :

    # read in analyzed InsetChart data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'All_Age_InsetChart.csv'))
    df['date'] = pd.to_datetime(df['date'])
    df = df.groupby(['date'] + sweep_variables)[channels_inset_chart].agg(np.mean).reset_index()

    # make InsetChart plot
    fig1 = plt.figure('InsetChart', figsize=(12,6))
    fig1.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig1.suptitle(f'Analyzer: InsetChartAnalyzer')
    axes = [fig1.add_subplot(2, 2, x + 1) for x in range(4)]
    for ch, channel in enumerate(channels_inset_chart) :
        ax = axes[ch]
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axes[-1].legend(title=sweep_variables)
    fig1.savefig(os.path.join(working_dir, expt_name, 'InsetChart.png'))


def plot_summary_report(sweep_variables, channels_summary_report=None, Uage='U5'):
    df = pd.read_csv(os.path.join(working_dir,expt_name, f'{Uage}_PfPR_ClinicalIncidence.csv'))
    df['date'] = df.apply(lambda x: datetime.date(int(x['year']), int(x['month']), 1), axis=1)
    df.columns = [x.replace(f' {Uage}', '') for x in df.columns]

    if channels_summary_report is None:
        channels_summary_report = ['Pop', 'Cases', 'Severe cases', 'PfPR']

    ## Aggregate runs and time (for simplicity take mean across all)!
    df = df.groupby(['date'] + sweep_variables)[channels_summary_report].agg(np.mean).reset_index()

    fig1 = plt.figure('MalariaSummaryReport', figsize=(12, 6))
    fig = plt.figure(figsize=(12, 6))
    fig.subplots_adjust(right=0.97, left=0.08, hspace=0.5, wspace=0.35, top=0.83, bottom=0.10)
    axes = [fig.add_subplot(2, 2, x + 1) for x in range(4)]
    fig.suptitle(f'MalariaSummaryReport {Uage}')

    for ai, channel in enumerate(channels_summary_report):
        ax = axes[ai]

        if channel == 'PfPR':
            ax.set_ylim(0, 1)
        else:
            ax.set_ylim(0, np.max(df[channel]))
        ax.set_xlabel('')

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axes[-1].legend(title=sweep_variables)
    fig.savefig(os.path.join(working_dir, expt_name, f'PfPR_ClinicalIncidence_{Uage}.png'))


def plot_events(event_list, sweep_variables) :

    # read in analyzed event data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'Event_Count.csv'))
    df['date'] = pd.to_datetime(df['date'])
    cov_channel_list = ['%s_Coverage' % x[9:] for x in event_list]
    cov_channel_list = [x for x in cov_channel_list if x in df.columns.values]
    df = df.groupby(['date'] + sweep_variables)[event_list + cov_channel_list].agg(np.mean).reset_index()

    # make event plot
    fig3 = plt.figure('Events', figsize=(12,3*len(event_list)))
    fig3.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig3.suptitle(f'Analyzer: ReceivedCampaignAnalyzer')
    axes = [fig3.add_subplot(len(event_list), 2, x + 1) for x in range(len(event_list))] #*2
    for ch, channel in enumerate(event_list) :
        ax = axes[ch] #*2
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

        coverage_channel = '%s_Coverage' % channel[9:]
        if coverage_channel not in df.columns.values :
            continue
        ax = axes[ch*2+1]
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[coverage_channel], '-', linewidth=0.8, label=p)
        ax.set_title(coverage_channel)
        ax.set_ylabel(coverage_channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    axes[-1].legend(title=sweep_variables)
    fig3.savefig(os.path.join(working_dir, expt_name, 'Events.png'))


if __name__ == "__main__":
    SetupParser.init()

    """Set sweep_variables and event_list as required depending on experiment"""
    sweep_variables = ['cm_cov_U5', 'itn_campaign', 'Run_Number']  #'itn_cooverage'
    event_list = ['Received_Treatment', 'Received_Severe_Treatment']
    event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']
    channels_inset_chart = ['Statistical Population', 'New Clinical Cases', 'Adult Vectors', 'Infected']

    analyzers = [InsetChartAnalyzer(expt_name=expt_name,
                                    working_dir=working_dir,
                                    channels=channels_inset_chart,
                                    sweep_variables=sweep_variables),
                 MonthlyPfPRAnalyzerU5(expt_name=expt_name,
                                          working_dir=working_dir,
                                          start_year=2020,
                                          sweep_variables=sweep_variables),
                 ReceivedCampaignAnalyzer(expt_name=expt_name,
                                          working_dir=working_dir,
                                          channels=event_list,
                                          start_year=2020,
                                          sweep_variables=sweep_variables),
                 ]

    am = AnalyzeManager(expt_id, analyzers=analyzers)
    am.analyze()

    sweep_vars_for_plotting = [x for x in sweep_variables if x != 'Run_Number']
    plot_inset_chart(channels_inset_chart, sweep_vars_for_plotting)
    plot_summary_report(sweep_vars_for_plotting)
    plot_events(event_list, sweep_vars_for_plotting)
    plt.show()