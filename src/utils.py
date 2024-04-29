"""
How Pairs of Insertion Mutations Impact Protein Structure: An Exhaustive Computational Study

2024

Changrui Li, Yang Zheng, and Filip Jagodzinski
"""

# Standard libary imports
import os

# Third party imports
import yaml, pickle

def load_config(config_file):
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

def inquiry(protein, indel1=None, aa1=None, indel2=None, aa2=None):
    """
    Search for mutations in a protein that match given criteria.
    If no specific criteria are provided, all mutations for the protein are returned.

    Args:
    protein (Protein): The protein object to search within.
    indel1 (str(int)): The position of the first insertion.
    aa1 (str): The amino acid of the first insertion.
    indel2 (str(int)): The position of the second insertion.
    aa2 (str): The amino acid of the second insertion.

    Returns:
    list of Mutation: A list of mutations that match the criteria.
    """
    results = []
    for mutation in protein.mutations:
        match = True
        if indel1 is not None and mutation.indel1 != indel1:
            match = False
        if aa1 is not None and mutation.aa1 != aa1:
            match = False
        if indel2 is not None and mutation.indel2 != indel2:
            match = False
        if aa2 is not None and mutation.aa2 != aa2:
            match = False
        if match:
            results.append(mutation)
    return results
