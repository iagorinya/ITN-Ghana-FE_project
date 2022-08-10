from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser

from analyzer_collection import *
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_example_w8b'
expt_id = '166b4fe3-e209-ed11-a9fa-b88303911bc1'  ## change expt_id
working_dir = os.path.join(os.getcwd(), 'simulation_outputs')


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
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axes[-1].legend(title=sweep_variables)
    fig1.savefig(os.path.join(working_dir, expt_name, 'InsetChart.png'))


def plot_summary_report(sweep_variables, ip_name, ipfilter):
    # read in analyzed summary reporrt
    channels_summary_report = ['Pop', 'Cases', 'Severe cases', 'PfPR']

    df = pd.DataFrame()
    for ipf in ipfilter:
        g = pd.read_csv(os.path.join(working_dir, expt_name, f'U5{ipf}_PfPR_ClinicalIncidence.csv'))
        g.columns = [col.replace(' U5', '') for col in g.columns]
        g[ip_name] = ipf.replace('_', '')
        df = pd.concat([g, df])

    df['date'] = df.apply(lambda x: datetime.date(int(x['year']), int(x['month']), 1), axis=1)

    # take mean over all Runs in report
    df = df.groupby(['date', 'month', 'year'] + sweep_variables + [ip_name])[channels_summary_report].agg(
        np.mean).reset_index()

    # make summary report plot
    fig2 = plt.figure('Summary Report', figsize=(6, 5))
    fig2.subplots_adjust(right=0.96, left=0.12, hspace=0.55, wspace=0.35, top=0.83, bottom=0.10)
    axes = [fig2.add_subplot(2, 2, x + 1) for x in range(4)]
    fig2.suptitle(f'Monthly_U5_PfPR with IP')

    for ai, channel in enumerate(channels_summary_report):
        ax = axes[ai]
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.set_ylim(0, max([1, 1.1 * np.max(df[channel])]))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

        for p, pdf in df.groupby(ip_name):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)

    axes[-1].legend(title=ip_name)
    fig2.savefig(os.path.join(working_dir, expt_name, 'U5_IP_PfPR_ClinicalIncidence.png'))


def plot_individual_events(event_list, ip_name, sweep_variables):
    # read in analyzed event data
    df = pd.read_csv(os.path.join(working_dir, expt_name, 'IndividualEvents_all_years.csv'))
    df = df.loc[df['Event_Name'].isin(event_list)]

    ## Filter to last year of simulation and one seed to reduce output
    # df = df.loc[df['Time'] >= max(df['Time']) - 365]
    df['date'] = df.apply(lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(x['Day'] - 1), axis=1)

    # 'date'
    adf = df.groupby(['Event_Name', 'Run_Number', ip_name] + sweep_variables)[['Individual_ID']].agg(
        'count').reset_index().rename(columns={'Individual_ID': 'count'})
    adf = adf.groupby(['Event_Name', ip_name] + sweep_variables)['count'].agg(np.mean).reset_index()
    adf.sort_values(by=['count'])

    # make event plot
    fig3 = plt.figure('Events', figsize=(12, 3 * len(event_list)))
    fig3.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig3.suptitle(f'Sum of individual events by {ip_name}')
    axes = [fig3.add_subplot(len(event_list), 2, x + 1) for x in range(len(event_list) * 2)]
    for ch, channel in enumerate(event_list):
        edf = adf.loc[adf['Event_Name'] == channel]
        ax = axes[ch]
        ax.set_title(channel)
        ax.set_ylabel(f'Count of {channel}')
        ax.set_xlabel(f'Access')
        for ipval in adf[ip_name].unique():
            pdf = edf.loc[edf[ip_name] == ipval]
            if not pdf.empty:
                ax.bar(pdf[ip_name].unique(), pdf['count'], alpha=0.4)

    fig3.savefig(os.path.join(working_dir, expt_name, 'Individual_Events.png'))


if __name__ == "__main__":
    # SetupParser.init()

    sweep_variables = ['itn_coverage', 'Run_Number']
    event_list = ['Received_ITN', 'Received_Treatment', 'Received_Severe_Treatment']
    channels_inset_chart = ['Statistical Population', 'New Clinical Cases', 'Adult Vectors', 'Infected']
    ip_name = 'Access'
    ipfilter = ['_accesslow', '_accesshigh']

    start_year = 2020
    # analyzers to run
    analyzers = [
        InsetChartAnalyzer(expt_name=expt_name,
                           working_dir=working_dir,
                           channels=channels_inset_chart,
                           sweep_variables=sweep_variables),
        MonthlyPfPRAnalyzerU5(expt_name=expt_name,
                              working_dir=working_dir,
                              start_year=start_year,
                              end_year=start_year + 1,
                              sweep_variables=sweep_variables),
        MonthlyPfPRAnalyzerU5IP(expt_name=expt_name,
                                working_dir=working_dir,
                                start_year=start_year,
                                end_year=start_year + 1,
                                sweep_variables=sweep_variables,
                                ipfilter='_accesslow'),
        MonthlyPfPRAnalyzerU5IP(expt_name=expt_name,
                                working_dir=working_dir,
                                start_year=start_year,
                                end_year=start_year + 1,
                                sweep_variables=sweep_variables,
                                ipfilter='_accesshigh'),
        IndividualEventsAnalyzer(expt_name=expt_name,
                                 working_dir=working_dir,
                                 start_year=start_year,
                                 sweep_variables=sweep_variables)
    ]
    am = AnalyzeManager(expt_id, analyzers=analyzers)
    am.analyze()

    sweep_vars_for_plotting = [x for x in sweep_variables if x != 'Run_Number']
    plot_inset_chart(channels_inset_chart, sweep_vars_for_plotting)
    plot_summary_report(sweep_vars_for_plotting, ip_name, ipfilter)
    plot_individual_events(event_list, ip_name, sweep_vars_for_plotting)
    plt.show()
