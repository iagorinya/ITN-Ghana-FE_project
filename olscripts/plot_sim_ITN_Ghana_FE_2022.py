import os
import datetime
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
palette = sns.color_palette("tab10")


def plot_inset_chart(channels_inset_chart, sweep_variables):
    # read in analyzed InsetChart data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'All_Age_InsetChart.csv'))
    df['date'] = pd.to_datetime(df['date'])

    ## Aggregate runs
    df = df.groupby(['date'] + sweep_variables)[channels_inset_chart].agg(np.mean).reset_index()

    # make InsetChart plot
    fig1 = plt.figure('InsetChart', figsize=(12, 6))
    fig1.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig1.suptitle(f'InsetChartAnalyzer')
    axes = [fig1.add_subplot(2, 2, x + 1) for x in range(4)]
    for ch, channel in enumerate(channels_inset_chart):
        ax = axes[ch]

        if channel == 'PfHRP2 Prevalence':
            ax.set_ylim(0, 1)
        else:
            ax.set_ylim(0, np.max(df[channel]))

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axes[-1].legend(title=sweep_variables)
    fig1.savefig(os.path.join(working_dir, expt_name, 'InsetChart.png'))


def plot_summary_report(sweep_variables, channels_summary_report=None, Uage='U5'):
    df = pd.read_csv(os.path.join(sim_dir, f'{Uage}_PfPR_ClinicalIncidence.csv'))
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


def plot_agebin_summary_report(sweep_variables, channels_summary_report=None):
    # read in analyzed summary reporrt
    if channels_summary_report is None:
        channels_summary_report = ['Pop', 'Cases', 'Severe cases', 'PfPR']
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'Agebin_PfPR_ClinicalIncidence_annual.csv'))
    df = df.sort_values(by='agebin')
    # take mean over all years in report
    df = df.groupby(['agebin'] + sweep_variables)[channels_summary_report].agg(np.mean).reset_index()

    # make summary report plot
    fig2 = plt.figure('Summary Report', figsize=(12, 6))
    fig2.subplots_adjust(right=0.96, left=0.12, hspace=0.55, wspace=0.35, top=0.83, bottom=0.10)
    axes = [fig2.add_subplot(2, 2, x + 1) for x in range(4)]
    fig2.suptitle(f'MalariaSummaryReport (agebin)')

    for ai, channel in enumerate(channels_summary_report):
        ax = axes[ai]
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.set_ylim(0, max([1, 1.1 * np.max(df[channel])]))
        ax.set_xlabel('Agebin')

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['agebin'], pdf[channel], '-', linewidth=0.8, label=p)

    axes[-1].legend(title=sweep_variables)
    fig2.savefig(os.path.join(working_dir, expt_name, 'Agebin_PfPR_ClinicalIncidence.png'))


def plot_events(event_list, sweep_variables):
    # read in analyzed event data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'Event_Count.csv'))
    df['date'] = pd.to_datetime(df['date'])
    cov_channel_list = [f'{x[9:]}_Coverage' for x in event_list]
    cov_channel_list = [x for x in cov_channel_list if x in df.columns.values]
    df = df.groupby(['date'] + sweep_variables)[event_list + cov_channel_list].agg(np.mean).reset_index()

    # make event plot
    fig3 = plt.figure('Events', figsize=(12, 3 * len(event_list)))
    fig3.subplots_adjust(right=0.97, left=0.08, hspace=0.5, wspace=0.35, top=0.83, bottom=0.10)
    fig3.suptitle(f'ReceivedCampaignAnalyzer')
    axes = [fig3.add_subplot(len(event_list), 2, x + 1) for x in range(len(event_list) * 2)]
    for ch, channel in enumerate(event_list):
        ax = axes[ch * 2]
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

        coverage_channel = f'{channel[9:]}_Coverage'
        if coverage_channel not in df.columns.values:
            continue
        ax = axes[ch * 2 + 1]
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[coverage_channel], '-', linewidth=0.8, label=p)
        ax.set_title(coverage_channel)
        ax.set_ylabel(coverage_channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    axes[-1].legend(title=sweep_variables)
    fig3.savefig(os.path.join(working_dir, expt_name, 'Events.png'))


def plot_transmission(sweep_variables, time_res='Monthly', selected_year=None):
    if selected_year is not None:
        selected_year = f'_{selected_year}'
    else:
        selected_year = '_all_years'

    df = pd.read_csv(os.path.join(sim_dir, f'{time_res}_transmission_report{selected_year}.csv'))
    x_var = 'Year'
    if 'date' in df.columns:
        df['date'] = pd.to_datetime(df['date'])
        x_var = 'date'

    channels = ['Daily Bites per Human', 'Daily EIR', 'Rainfall', 'PfHRP2 Prevalence']
    channels = [x.replace('Daily', time_res) for x in channels]

    df = df.groupby([x_var] + sweep_variables)[channels].agg(np.mean).reset_index()

    fig = plt.figure(figsize=(12, 6))
    fig.subplots_adjust(right=0.97, left=0.08, hspace=0.5, wspace=0.35, top=0.83, bottom=0.10)
    axes = [fig.add_subplot(2, 2, x + 1) for x in range(4)]
    fig.suptitle(f'TransmissionReport')

    for ai, channel in enumerate(channels):
        ax = axes[ai]
        ax.set_title(channel)
        ax.set_ylabel(channel)
        if channel == 'PfHRP2 Prevalence':
            ax.set_ylim(0, 1)
        else:
            ax.set_ylim(0, np.max(df[channel]))
        ax.set_xlabel('Time')
        if x_var == 'date':
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m\n%y'))

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf[x_var], pdf[channel], '-', linewidth=0.8, label=p)

    axes[-1].legend(title=sweep_variables)
    fig.savefig(os.path.join(working_dir, expt_name, f'Transmission_{time_res}.png'))


if __name__ == "__main__":
    user = os.getlogin()  # user initials
    expt_name = f'{user}_FE_2022_example_w4'
    working_dir = os.path.join('simulation_outputs')
    sim_dir = os.path.join(working_dir, expt_name)

    """Set sweep_variables and event_list as required depending on experiment"""
    sweep_variables = ['cm_cov_U5', 'itn_coverage', 'smc_coverage']
    event_list = ['Received_Treatment', 'Received_ITN', 'Received_SMC']
    channels_inset_chart = ['Statistical Population', 'New Clinical Cases', 'Adult Vectors', 'Infected']

    """Generate plots"""
    plot_inset_chart(channels_inset_chart, sweep_variables)
    plot_summary_report(sweep_variables, Uage='U5')  # (MalariaSummaryReport aggregated age)
    plot_agebin_summary_report(sweep_variables)  # (MalariaSummaryReport agebins)
    plot_events(event_list, sweep_variables)
    # plot_transmission(sweep_variables) # optional
