"""
Identifying And Comparison Analyze Impactful Pairs of Insertion Mutations in Proteins

2024

Changrui Li and Yang Zheng
"""

# standard library imports
import os


# import - third paty
import json, pickle, statistics,math

# import - local
from utils import load_config, get_filenames, load_protein, crop_mutations
from protein import Protein, Mutation


def main():
    # main yaml path
    main_yaml_path = '../config/main.yaml'
    
    # settings configurations
    config = load_config(main_yaml_path)
    metrics = config.get('metrics')
    proteins = config.get('proteins')
    saved_protein_data_path = config.get('saved_protein_data_path')
    
    # read filenames
    filenames = list(get_filenames(saved_protein_data_path))
    
    # exhuastive data for three proteins
    proteins_data = {filename[:-4]:load_protein(os.path.join(saved_protein_data_path,filename)) for filename in filenames}
    print("# protein mutations: exhuastive set")
    for protein_name, data in proteins_data.items():
        print(data)
    print()


    # crop mutations for mutations for three proteins
    cropped_mutations_groups = {}
    for protein in proteins:
        cropped_mutations_groups[protein] = crop_mutations(proteins_data[protein].mutations)

    # metrics
    print("# metrics")
    [print(metric) for metric in metrics]
    print()

    # cropped mutations
    print("# cropped mutations")
    for protein_name, mutations in cropped_mutations_groups.items():
        print(f"<protein name={protein_name}, number of cropped mutations={len(mutations)}>")
    print()

    # mutation sample
    print("# mutation sample")
    print(cropped_mutations_groups['1hhp'][0])

if __name__ == "__main__":
    main()
    pass