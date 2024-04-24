"""
Identifying And Comparison Analyze Impactful Pairs of Insertion Mutations in Proteins

2024

Changrui Li and Yang Zheng
"""

# standard library imports
import os, sys


# import - third paty
import json, pickle, statistics,math,statistics
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
# import wandb

# import - local
from utils import load_config, get_filenames,load_protein, crop_mutations
from protein import Protein, Mutation, _mutants_outside_two_stdev_for_metric
from data.amino_acids import sorted_amino_acid_list
# discover the average hbond that grouped by amino acid type

# import tasks
from amino_acid_groups_over_three_metrics_table import amino_acid_groups_over_three_metrics_table

def tasks(cropped_mutations_groups):
    # amino_acid_groups_over_three_metrics_table()
    
    # two standard deviation analyze, there are multiple tasks within below function, go there and check for more detail
    from outlier_two_standard_deviation import outlier_two_stand_deviation
    outlier_two_stand_deviation(cropped_mutations_groups)

def main():
    # main yaml path
    main_yaml_path = '../config/main.yaml'
    
    # settings and retrive the saved_protein_data_path
    config = load_config(main_yaml_path)
    metrics = config.get('metrics')
    saved_protein_data_path = config.get('saved_protein_data_path')
    
    # read filenames
    filenames = list(get_filenames(saved_protein_data_path))
    
    # set protein data
    proteins_data = {filename[:-4]:load_protein(os.path.join(saved_protein_data_path,filename)) for filename in filenames}

    # crop mutations for mutations of each protein
    cropped_mutations_groups = {} 
    cropped_mutations_groups['1hhp'] = crop_mutations(proteins_data['1hhp'].mutations)
    cropped_mutations_groups['1csp'] = crop_mutations(proteins_data['1csp'].mutations)
    cropped_mutations_groups['1crn'] = crop_mutations(proteins_data['1crn'].mutations)

    # end of data prepare


    # running tasks
    tasks(cropped_mutations_groups)
    


    


    
if __name__ == "__main__":
    main()
    pass