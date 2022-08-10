from analyzer_collection import *
import os
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_w7pick_up50'  ## change w6a or w6b
expt_id = '22298e88-2018-ed11-a9fb-b88303911bc1'  ## change expt_id
working_dir = os.path.join('simulation_outputs')

serialize_years = 50  # Same as in run_exampleSim_w6a.py
step = 'pickup'  # 'pickup'  # 'burnin'


def plot_inset_chart(channels_inset_chart, sweep_variables):
    # read in analyzed InsetChart data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'All_Age_InsetChart.csv'))
    df['date'] = pd.to_datetime(df['date'])
    df = df.groupby(['date'] + sweep_variables)[channels_inset_chart].agg(np.mean).reset_index()

    # make InsetChart plot
    fig1 = plt.figure('InsetChart', figsize=(12, 6))
    fig1.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig1.suptitle(f'Analyzer: InsetChartAnalyzer')
    axes = [fig1.add_subplot(2, 2, x + 1) for x in range(4)]
    for ch, channel in enumerate(channels_inset_chart):
        ax = axes[ch]
        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    if len(sweep_variables) > 1:
        axes[-1].legend(title=sweep_variables)
    fig1.savefig(os.path.join(working_dir, expt_name, 'InsetChart.png'))


def plot_summary_report(sweep_variables, selected_age=5):
    # read in analyzed summary reporrt
    channels_summary_report = ['Pop', 'Cases', 'Severe cases', 'PfPR']
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'Agebin_PfPR_ClinicalIncidence_annual.csv'))
    df = df.sort_values(by='agebin')
    df = df.loc[df['agebin'] == selected_age]  # select single agebin
    ## Take mean of runs
    df = df.groupby(['year'] + sweep_variables)[channels_summary_report].agg(np.mean).reset_index()

    # make summary report plot
    fig2 = plt.figure('Summary Report', figsize=(6, 5))
    fig2.subplots_adjust(right=0.96, left=0.12, hspace=0.55, wspace=0.35, top=0.83, bottom=0.10)
    axes = [fig2.add_subplot(2, 2, x + 1) for x in range(4)]
    fig2.suptitle(f'Pick up simulation - annual summary report')

    for ai, channel in enumerate(channels_summary_report):
        ax = axes[ai]
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.set_ylim(0, max([1, np.max(df[channel])]))
        ax.set_xlabel('year')

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['year'], pdf[channel], '-', linewidth=0.8, label=p)

    axes[-1].legend(title=sweep_variables)
    fig2.savefig(os.path.join(working_dir, expt_name, 'PfPR_ClinicalIncidence.png'))


if __name__ == "__main__":
    SetupParser.init()

    """Set sweep_variables and event_list as required depending on experiment"""
    channels_inset_chart = ['Statistical Population', 'New Clinical Cases', 'Adult Vectors', 'Infected']
    sweep_variables = ['Run_Number']
    if step == 'pickup':
        sweep_variables = ['Run_Number', 'itn_coverage']
        # sweep_variables = ['Run_Number', 'cm_cov_U5', 'itn_coverage']

    analyzers_burnin = [InsetChartAnalyzer(expt_name=expt_name,
                                           working_dir=working_dir,
                                           channels=channels_inset_chart,
                                           start_year=1980 - serialize_years,
                                           sweep_variables=sweep_variables),
                        ]

    analyzers_pickup = [InsetChartAnalyzer(expt_name=expt_name,
                                           working_dir=working_dir,
                                           channels=channels_inset_chart,
                                           start_year=1980,
                                           sweep_variables=sweep_variables),
                        # AnnualAgebinPfPRAnalyzer(expt_name=expt_name,
                        #                          working_dir=working_dir,
                        #                          start_year=1961,
                        #                          sweep_variables=sweep_variables)
                        ]

    if step == 'burnin':
        am = AnalyzeManager(expt_id, analyzers=analyzers_burnin)
        #am.analyze()
        sweep_vars_for_plotting = ['Run_Number']
        plot_inset_chart(channels_inset_chart, sweep_vars_for_plotting)

    elif step == 'pickup':
        am = AnalyzeManager(expt_id, analyzers=analyzers_pickup)
        am.analyze()
        sweep_vars_for_plotting = [x for x in sweep_variables if x != 'Run_Number']
        plot_inset_chart(channels_inset_chart, sweep_vars_for_plotting)
        #plot_summary_report(sweep_vars_for_plotting)
    else:
        print('Please define step, options are burnin or pickup')

plt.show()
