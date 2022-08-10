import os
import datetime
import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser

SetupParser.default_block = 'HPC'

user = os.getlogin()  # user initials
expt_name = f'{user}_FE_2022_burnin10'
expt_id = '5784135b-2214-ed11-a9fb-b88303911bc1'  ## change expt_id
working_dir = os.path.join('simulation_outputs')


class VectorSpeciesAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, sweep_variables=None, channels=None, working_dir=".", start_year=2022):
        super(VectorSpeciesAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/VectorSpeciesReport.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.inset_channels = channels or ['Adult Vectors Per Node']
        self.expt_name = expt_name
        self.start_year = start_year

    def select_simulation_data(self, data, simulation):
        species = data[self.filenames[0]]['Header']['Subchannel_Metadata']['MeaningPerAxis'][0][0]
        simdata = pd.DataFrame({'%s %s' % (x, sp): data[self.filenames[0]]['Channels'][x]['Data'][s] for x in self.inset_channels for s, sp in enumerate(species)})
        simdata['Time'] = simdata.index
        simdata['Day'] = simdata['Time'] % 365
        simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(int(x['Day']) - 1), axis=1)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            elif sweep_var == 'Run_Number' :
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        species = [x.split(' ')[-1] for x in adf.columns.values if self.inset_channels[0] in x]

        adf['total vectors'] = adf['%s %s' % (self.inset_channels[0], species[0])]
        for i in range(1, len(species)) :
            adf['total vectors'] += adf['%s %s' % (self.inset_channels[0], species[i])]
        for sp in species :
            adf['normalized %s' % sp] = adf['%s %s' % (self.inset_channels[0], sp)]/adf['total vectors']

        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'Vector_Species.csv'), index=False)


if __name__ == "__main__":
    SetupParser.init()

    # set desired channels to analyze and plot
    channels = ['Adult Vectors Per Node']

    # call analyzer to grab EMOD output and process into analyzed data
    analyzers = [VectorSpeciesAnalyzer(expt_name=expt_name,
                                       channels=channels,
                                       working_dir=working_dir)]

    am = AnalyzeManager(expt_id, analyzers=analyzers)
    am.analyze()