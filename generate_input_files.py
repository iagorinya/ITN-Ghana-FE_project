import os
import pandas as pd
import json
from dtk.tools.demographics.DemographicsGeneratorConcern import WorldBankBirthRateConcern, \
    EquilibriumAgeDistributionConcern, DefaultIndividualAttributesConcern
from dtk.tools.demographics.DemographicsGenerator import DemographicsGenerator
from dtk.tools.climate.ClimateGenerator import ClimateGenerator


def generate_demographics(demo_df, my_ds, demo_fname):
    demo_df = demo_df.loc[demo_df['repDS'] == my_ds]
    demo_df['population'] = 1000  # demo_df['population']/100

    br_concern = WorldBankBirthRateConcern(country="Ghana", birthrate_year=2016)

    chain = [
        DefaultIndividualAttributesConcern(),
        br_concern,
        EquilibriumAgeDistributionConcern(default_birth_rate=br_concern.default_birth_rate),
    ]

    current = DemographicsGenerator.from_dataframe(demo_df,
                                                   population_column_name='population',
                                                   latitude_column_name='lat',
                                                   longitude_column_name='lon',
                                                   node_id_from_lat_long=False,
                                                   concerns=chain,
                                                   load_other_columns_as_attributes=True,
                                                   include_columns=['repDS']
                                                   )
    current['Nodes'][0]['NodeID'] = 1
    with open(demo_fname, 'w') as fout:
        json.dump(current, fout, sort_keys=True, indent=4, separators=(',', ': '))


def generate_climate(demo_fname, my_ds):
    from simtools.SetupParser import SetupParser

    if not SetupParser.initialized:
        SetupParser.init('HPC')

    cg = ClimateGenerator(demographics_file_path=demo_fname, work_order_path='./wo.json',
                          climate_files_output_path=os.path.join(inputs_path, my_ds),
                          climate_project='IDM-Ghana',
                          start_year='2001', num_years='16')
    cg.generate_climate_files()


if __name__ == '__main__':
    inputs_path = os.path.join('input/Ghana')
    if not os.path.exists(inputs_path):
        os.mkdir(inputs_path)

    ##FIXME read in csv with ecozone names and centroids for long/lat
    df = pd.DataFrame(data={'repDS': ['GHA_ecozone1', 'GHA_ecozone2', 'GHA_ecozone3'],
                            'population': [1000, 1000, 1000], ## dummy data
                            'lat': [5.760759295941768, 5.760759295941768, 5.760759295941768],  ## dummy data
                            'lon': [-0.4473415119456551, -0.4473415119456551, -0.4473415119456551]})  ## dummy data

    repDS_list = list(df['repDS'].unique())
    for my_ds in repDS_list:
        print(my_ds)
        if not os.path.exists(os.path.join(inputs_path, my_ds)):
            os.makedirs(os.path.join(inputs_path, my_ds))

        demo_fname = os.path.join(inputs_path, my_ds, f'{my_ds}_demographics.json')
        generate_demographics(df, my_ds, demo_fname)
        generate_climate(demo_fname, my_ds)
