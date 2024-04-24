import os, json, pickle, statistics
from multiprocessing import Pool


def _extract_metric(mutations, metric_name):
    """Helper method to extract values of a metric across all mutations."""
    return [getattr(mutation, metric_name) for mutation in mutations]

def _mutants_outside_two_stdev_for_metric(mutations, metric_name, get_distance=None):
    """Helper method to get mutations outside two standard deviations for a specific metric."""
    """this distance means that if calculating how far from the mean """
    values = _extract_metric(mutations, metric_name)
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)

    if get_distance==True:
        
        local_mutations = []
        for mutation in mutations:
            value = abs(getattr(mutation, metric_name) - mean)
            if value > 2 * stdev:
                setattr(mutation,metric_name,value)
            local_mutations.append(mutation)
        return local_mutations
    else:
        return [mutation for mutation in mutations if abs(getattr(mutation, metric_name) - mean) > 2 * stdev]

class Protein:
    def __init__(self, name):
        self.name = name
        self.mutations = []
        self.helix = []
        self.sheet = []
    def add_mutation(self, mutation):
        """Add a mutation to the protein's list of mutations."""
        self.mutations.append(mutation)

    def group_mutations_by_amino_acid(self):
        """Group mutations by hydrogen bond counts."""
        groups = {}
        for mutation in self.mutations:
            aa1 = mutation.aa1
            aa2 = mutation.aa2
            if aa1 not in groups:
                groups[aa1] = []
            if aa2 not in groups:
                groups[aa2] = []
            groups[aa1].append(mutation)
            groups[aa2].append(mutation)
        return groups

class Mutation:
    def __init__(self, indel1, aa1, indel2, aa2, hbond_count, size_of_largest_cluster, rigidity_order_parameter, cluster_configuration_entropy, hbond_diff=None):
        self.indel1 = indel1
        self.aa1 = aa1
        self.indel2 = indel2
        self.aa2 = aa2
        self.hbond_count = hbond_count
        self.size_of_largest_cluster = size_of_largest_cluster
        self.rigidity_order_parameter = rigidity_order_parameter
        self.cluster_configuration_entropy = cluster_configuration_entropy
        self.hbond_diff = hbond_diff

    def __repr__(self):
        return f"<Mutation indel1={self.indel1}, aa1={self.aa1}, indel2={self.indel2}, aa2={self.aa2}, hbond={self.hbond_count}, size of largest cluster={self.size_of_largest_cluster}, rigidity order parameter={self.rigidity_order_parameter}, cluster configuration entropy={self.cluster_configuration_entropy}>"


# class definitions for Protein and Mutation here ...

def save_protein_to_file(protein_data, local_path, filename):
    filename = os.path.join(local_path,filename)
    with open(filename, 'wb') as file:
        pickle.dump(protein_data, file)

# load protein class
def load_protein_from_file(filename):
    with open(filename, 'rb') as file:
        return pickle.load(file)

def load_mutation_from_file(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
        
        # Extracting relevant parameters from the JSON data
        protein_name = data["params0"]
        action = data["params1"]
        indel1 = data["params2"]
        aa1 = data["params3"]
        indel2 = data["params4"]
        aa2 = data["params5"]
        hbond_count = data["hbond_count"]
        # Extracting metrics from the nested 'kinari_metrics' key
        metrics = data["kinari_metrics"]
        size_of_largest_cluster = metrics["size of largest Clust"]
        rigidity_order_parameter = metrics["Rigidity order parameter"]
        cluster_configuration_entropy = metrics["Cluster configuration entropy"]
        # i forget the number of hbond
        
        return Mutation(indel1, aa1, indel2, aa2, hbond_count, size_of_largest_cluster, rigidity_order_parameter, cluster_configuration_entropy)


def load_mutations_from_subfolder(subfolder_path):
    mutations = []

    # Check if subfolder contains files
    if not os.listdir(subfolder_path):
        return mutations  # Return empty list if subfolder has no files
    
    # Otherwise, iterate over each file in the subfolder
    for filename in os.listdir(subfolder_path):
        full_path = os.path.join(subfolder_path, filename)
        if os.path.isfile(full_path):
            mutation = load_mutation_from_file(full_path)
            mutations.append(mutation)
    
    return mutations

def load_protein_from_directory(protein_action_folder):
    protein_name = os.path.basename(os.path.dirname(protein_action_folder))
    protein = Protein(protein_name)

    target_depth = len(protein_action_folder.split(os.path.sep)) + 2

    for root, _, files in os.walk(protein_action_folder):
        current_depth = len(root.split(os.path.sep))
        if current_depth > target_depth:
            break  # Don't traverse deeper than necessary

        if current_depth == target_depth:
            with Pool(processes=28) as pool:
                mutations = pool.map(load_mutation_from_file, [os.path.join(root, file) for file in files])
                for mutation in mutations:
                    protein.add_mutation(mutation)
                    print(f'Added: {mutation}')

    return protein


def main(protein_name):
    action = 'ins'
    root_path = os.path.join('/research/jagodzinski/lic7/data/mutation_results',protein_name,action)
    print(f'Start loading: {protein_name}.')
    protein_class_data = load_protein_from_directory(root_path)
    print(f'End loading: {protein_name}.')
    
    print(f'Start saving: {protein_name}.pkl .')
    protein_class_path = '/research/jagodzinski/lic7/data/save_classes/saved_with_all_patterns'
    # os.makedirs(protein_class_path,exist_ok=True)
    save_protein_to_file(protein_class_data,protein_class_path,protein_name+'.pkl')
    print(f'End saving: {protein_name}.pkl .')

def load(protein_name):
    protein_class_path = '/research/jagodzinski/lic7/data/save_classes/saved_with_all_patterns'
    # load protein class
    protein_class_data = load_protein_from_file(os.path.join(protein_class_path,protein_name+'.pkl'))
    # print(len(protein_class_data.mutations))
    # print(type(protein_class_data))
    mutations = protein_class_data.mutations
    
    three_metrix = ['hbond_count','size_of_largest_cluster','rigidity_order_parameter','cluster_configuration_entropy']
    for metrix in three_metrix:
        print( metrix, len(_mutants_outside_two_stdev_for_metric(mutations, metrix)))


 
if __name__ == '__main__':
    import sys
    action = sys.argv[1]
    protein = sys.argv[2]
    if action == 'g':
        main(protein)
    if action == 'l':
        load(protein)
    pass

