import matplotlib, os, sys, wandb, math, statistics, requests
from statistics import mean, stdev
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import numpy as np
from collections import Counter, defaultdict
from protein_class import load_protein_from_file, _mutants_outside_two_stdev_for_metric, Protein, Mutation
from multiprocessing import Pool
from main import get_sorted_two_AA_list_hydrophobicity, get_sorted_two_AA_lines, create_folder, upload_to_wandb
from scipy.stats import norm, ttest_ind
from itertools import combinations
plt.rcParams.update({'font.size': 14 * 1.8})

def main_load(protein_name):
    protein_class_path = '/research/jagodzinski/lic7/data/save_classes/saved_with_all_patterns'
    pdb_url = f'https://files.rcsb.org/download/{protein_name}.pdb'
    response = requests.get(pdb_url)
    pdb_content = response.text

    global plot_type
    plot_type = 'scatter'
    # load mutation class
    protein = load_protein_from_file(os.path.join(protein_class_path, protein_name + '.pkl'))
    mutations = protein.mutations
    cropped_mutations = [mutation for mutation in mutations if abs(int(mutation.indel1) - int(mutation.indel2)) > 12]

    helix, sheet = parse_pdb(protein, pdb_content)
    protein.helix = helix
    protein.sheet = sheet

    # check if d1 and d2 - 1 is in alpha helix sequence
    upload = False
    grouped_mutations = group_by_insertion_position(helix, cropped_mutations)
    metrics_list = ['hbond_count', 'size_of_largest_cluster', 'rigidity_order_parameter', 'cluster_configuration_entropy']
    compare_significance(grouped_mutations, metrics_list, protein_name, upload)

    exit()

# parse pdb file and return the helix and sheet sequence
def parse_pdb(protein, pdb_content):
    helix = []
    sheet = []
    for line in pdb_content.split('\n'):
        # only need chain A
        if 'HELIX' in line and line.split()[4] == 'A':
           line = line.split()
           helix.append((int(line[5]), int(line[8])))
        if 'SHEET' in line and line.split()[2] == 'A':
           line = line.split()
           sheet.append((int(line[6]), int(line[9])))
    helix = sorted(helix, key=lambda x: int(x[0]))
    sheet = sorted(sheet, key=lambda x: int(x[0]))
    return helix, sheet

# take a list of mutations, group the mutations by the insertion position, [both in sequence], [one in sequence], [none in sequence]
def group_by_insertion_position(sequence, mutations):
    grouped_mutations = [[], [], []]
    for mutant in mutations:
        index = 0
        for seq in sequence:
            # because after the first insertion, the second insertion will be +1, so we need to -1 to get the correct position
            if int(mutant.indel1) in range(seq[0], seq[1] + 1):
                index += 1
            if int(mutant.indel2) - 1 in range(seq[0], seq[1] + 1):
                index += 1
            if index == 2:
                break
        grouped_mutations[index].append(mutant)
    
    # print(len(grouped_mutations[2]))
    print(len(grouped_mutations[0]), len(grouped_mutations[1]), len(grouped_mutations[2]))
    return grouped_mutations
    # [print((mutation.indel1, mutation.indel2)) for mutation in grouped_mutations[2]]

# two sample t-test for each metric for group NN, 1N
def compare_significance(grouped_mutations, metrics_list, protein_name, upload):

    for metric in metrics_list:
        H_00 = [getattr(mutation, metric) for mutation in grouped_mutations[0]]
        H_01 = [getattr(mutation, metric) for mutation in grouped_mutations[1]]
        H_11 = [getattr(mutation, metric) for mutation in grouped_mutations[2]]
        
        groups = {0:H_00, 1:H_01, 2:H_11}
        group_names = [r'$I_{XX}$', r'$I_{XH}$', r'$I_{HH}$']
        permutations = [0, 1]
        if protein_name == '1crn':
            permutations.append(2)

        permutation_groups = list(combinations(permutations, 2))

        for group in permutation_groups: # eg: [0, 1]
            print(f'group: {group_names[group[0]]}, {group_names[group[1]]}')
            # print(group)
            group_1 = groups[group[0]]
            group_2 = groups[group[1]]
            # Perform two-sample t-test
            t_stat, p_val = ttest_ind(group_1, group_2)
            # calculate cohens_d
            cohens_d = (mean(group_1) - mean(group_2)) / (math.sqrt((stdev(group_1) ** 2 + stdev(group_2) ** 2) / 2))

            # Fit a normal distribution to group_1 and group_2 data
            mu1, std1 = norm.fit(group_1)
            mu2, std2 = norm.fit(group_2)

            plt.figure(figsize=(15, 10))  # Adjust the figure size as needed

            # Set up plot range
            xmin = min(mu1 - 3*std1, mu2 - 3*std2)
            xmax = max(mu1 + 3*std1, mu2 + 3*std2)
            x = np.linspace(xmin, xmax, 100)
            plt.xlim(xmin, xmax * 1.5)

            # Plot the PDF for Group group_1
            p1 = norm.pdf(x, mu1, std1)
            plt.plot(x, p1, 'blue', linewidth=2, label=f'Group {group_names[group[0]]} (mean={mu1:.3f})')
            plt.axvline(mu1, color='blue', linestyle='dashed', linewidth=2)
            plt.axvline(mu1 - std1, color='blue', linestyle='-.', linewidth=2)
            plt.axvline(mu1 + std1, color='blue', linestyle='-.', linewidth=2)

            # Plot the PDF for Group group_2
            p2 = norm.pdf(x, mu2, std2)
            plt.plot(x, p2, 'green', linewidth=2, label=f'Group {group_names[group[1]]} (mean={mu2:.3f})')
            plt.plot([], [], ' ', label=f'cohens d = {cohens_d:.3f}')
            plt.plot([], [], ' ', label=f'{group_names[group[0]]} sample size= {len(group_1)}')
            plt.plot([], [], ' ', label=f'{group_names[group[1]]} sample size= {len(group_2)}')
            
            plt.axvline(mu2, color='green', linestyle='dashed', linewidth=2)
            plt.axvline(mu2 - std2, color='green', linestyle='-.', linewidth=2)
            plt.axvline(mu2 + std2, color='green', linestyle='-.', linewidth=2)

            # Plot formatting
            title = f'{protein_name}_two-sample_t-test_{metric}_{group_names[group[0]]}_{group_names[group[1]]}'
            # plt.title({title})
            plt.xlabel(f'{metric}')
            plt.ylabel('Density')
            plt.legend(loc='upper right')

            ax = plt.gca()  # Gets the current axes.
            # Use the axes transforms to place a box around the actual plot area
            trans = ax.transAxes  # This is the key to use axes-relative coordinates

            # Create a new rectangle that exactly overlays the axes' position
            rect = Rectangle(
                (0, 0), 1, 1, transform=trans, 
                linewidth=3,  # Adjust the thickness of the border here
                edgecolor='black', facecolor='none', zorder=10
            )

            # Add the rectangle to the axes so it is drawn in axes-relative coordinates
            ax.add_patch(rect)

            # Show the plot
            path = f'p-test/{protein_name}_{metric}_{group_names[group[0]]}_{group_names[group[1]]}'
            plt.tight_layout()
            plt.savefig(path, dpi=400)
            if upload:
                upload_to_wandb(plt, f't-test of group {group_names[group[0]]}, {group_names[group[1]]}', 't-test per metric', title)

def run_program(protein_name, metric, value, operator):
    # Placeholder for your program's logic
    print(f"Running program with protein name: {protein_name}, metric: {metric}, value: {value}, operator: {operator}")

if __name__ == "__main__":
    try:
        protein_name = sys.argv[1]
    except:
        print('python3 create_plots.py 1hhp 0 4 =')
        print(['hbond_count','size_of_largest_cluster','rigidity_order_parameter','cluster_configuration_entropy'])
        exit()
    main_load(protein_name)
    


