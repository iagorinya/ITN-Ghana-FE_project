import os
import pandas as pd
import json
from dtk.tools.demographics.DemographicsGeneratorConcern import WorldBankBirthRateConcern, \
    EquilibriumAgeDistributionConcern, DefaultIndividualAttributesConcern
from dtk.tools.demographics.DemographicsGenerator import DemographicsGenerator
from dtk.tools.climate.ClimateGenerator import ClimateGenerator


def generate_demographics(demo_df, demo_fname):
    # Get WorldBank birth rate estimate
    # UPDATE country and birthrate_year if needed
    br_concern = WorldBankBirthRateConcern(country="Ghana", birthrate_year=2016)

    chain = [
        DefaultIndividualAttributesConcern(),
        br_concern,
        EquilibriumAgeDistributionConcern(default_birth_rate=br_concern.default_birth_rate),
    ]

    current = DemographicsGenerator.from_dataframe(demo_df,
                                                   population_column_name='population',
                                                   nodeid_column_name='nodeid',
                                                   node_id_from_lat_long=False,
                                                   concerns=chain,
                                                   load_other_columns_as_attributes=True,
                                                   include_columns=['Ecozone', 'rank'])  # Add any "optional" columns

    with open(demo_fname, 'w') as fout:
        json.dump(current, fout, sort_keys=True, indent=4, separators=(',', ': '))


def generate_climate(demo_fname):
    from simtools.SetupParser import SetupParser

    if not SetupParser.initialized:
        SetupParser.init('HPC')

    cg = ClimateGenerator(demographics_file_path=demo_fname, work_order_path='./wo.json',
                          climate_files_output_path=inputs_path,
                          climate_project='IDM-Ghana',
                          start_year='2001', num_years='16')
    cg.generate_climate_files()


if __name__ == '__main__':
    inputs_path = os.path.join('./', 'input/Ghana')
    if not os.path.exists(inputs_path):
        os.mkdir(inputs_path)

    df = pd.DataFrame(data={'nodeid': [1], 'population': [1000], 'ecozone': ['forest zone'],
                            'lat': [6.6666], 'lon': [-1.616271]})

    demo_fname = os.path.join(inputs_path, 'Ghana_demographics_Zone_2.json')
    generate_demographics(df, demo_fname)
    generate_climate(demo_fname)
