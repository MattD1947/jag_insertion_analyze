"""
Identifying And Comparison Analyze Impactful Pairs of Insertion Mutations in Proteins

2024

Changrui Li and Yang Zheng
"""

# Standard libary imports
import os

# Third party imports
import yaml, json, pickle, statistics

def load_config(config_file):
    """
        Reads the config yaml file and returns its contents.
    """
    with open(config_file,'r') as f:
        config = yaml.safe_load(f)
    f.close()

    return config

# load protein class
def load_protein(filename):
    with open(filename,'rb') as file:
        return pickle.load(file)

# get filenames
def get_filenames(target_path):
    for filename in os.listdir(path=target_path):
        yield filename 

#crop the mutations that two inserted positions that are closer than or equal to 12 since at 12 there is majority missing mutations that can not be generated 
def crop_mutations(mutations):
    return [mutation for mutation in mutations if abs(int(mutation.indel1) - int(mutation.indel2)) > 12]