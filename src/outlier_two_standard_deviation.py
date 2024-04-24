import pandas as pd

def task_group_mutations_by_amino_acid(proteins,metrics,filtered_mutations):
    # group mutations by amino acid
    from utils import get_sorted_amino_acid_list, single_frequency_bar_plot

    upload = False
    plot_key = project_name = 'outliers_two_standard_deviation_amino_acid_group_frequency' 
    x_label = 'amino acid'
    y_label = 'frequency'
    
    sorted_amino_acid_dict = {item[0]:0 for item in get_sorted_amino_acid_list()}
    
    count_mutations_by_amino_acid_group = {protein:{metric:sorted_amino_acid_dict.copy() for metric in metrics} for protein in proteins}
    for protein in proteins:
        for metric, mutations in filtered_mutations[protein].items():
            for mutation in mutations:
                # update inserted animo acid count
                count_mutations_by_amino_acid_group[protein][metric][mutation.aa1] += 1
                count_mutations_by_amino_acid_group[protein][metric][mutation.aa2] += 1

            x = [item[0] for item in get_sorted_amino_acid_list()]
            y = [item[1] for item in count_mutations_by_amino_acid_group[protein][metric].items()]
            top_5 = sorted(count_mutations_by_amino_acid_group[protein][metric].items(),key=lambda x:x[1],reverse=True)[:5]
            
            # set dummy label
            dummy_label = ''
            for item in top_5:
                dummy_label += " ".join([item[0]+":",str(item[1])+"\n"])
            dummy_label = dummy_label[:-1]  

            # Create a dummy handle
            from matplotlib.lines import Line2D
            dummy_handle = Line2D([], [], color='none', label=dummy_label)
            
            plot_title = "_".join([protein,metric,plot_key])

            single_frequency_bar_plot(x,y,x_label,y_label,dummy_handle,plot_title,upload,project_name,plot_key)

def task_group_mutations_by_position(proteins,metrics,filtered_mutations):
    # group mutations by position
    from utils import get_protein_chain_length, single_frequency_bar_plot

    upload = False
    plot_key = project_name = 'outliers_two_standard_deviation_position_group_frequency' 
    x_label = 'position'
    y_label = 'frequency'
    count_mutations_by_position_group = {protein:{metric:{} for metric in metrics} for protein in proteins}
    for protein in proteins:
        protein_chain_length = get_protein_chain_length(protein)
        sorted_position_dict = { str(i):0 for i in range(1,protein_chain_length+1+2)}
        for metric, mutations in filtered_mutations[protein].items():
            count_mutations_by_position_group[protein][metric] = sorted_position_dict
            for mutation in mutations:
                # update inserted position count
                count_mutations_by_position_group[protein][metric][mutation.indel1] += 1
                count_mutations_by_position_group[protein][metric][mutation.indel2] += 1
            
            x = [str(position) for position in range(1,protein_chain_length+1+2)]
            y = [item[1] for item in count_mutations_by_position_group[protein][metric].items()]
            top_5 = sorted(count_mutations_by_position_group[protein][metric].items(),key=lambda x:x[1],reverse=True)[:5]
            
            # set dummy label
            dummy_label = ''
            for item in top_5:
                dummy_label += " ".join([item[0]+":",str(item[1])+"\n"])
            dummy_label = dummy_label[:-1]  

            # Create a dummy handle
            from matplotlib.lines import Line2D
            dummy_handle = Line2D([], [], color='none', label=dummy_label)
            
            plot_title = "_".join([protein,metric,plot_key])

            single_frequency_bar_plot(x,y,x_label,y_label,dummy_handle,plot_title,upload,project_name,plot_key)

def task_group_mutations_by_amino_acid_pairs(proteins,metrics,filtered_mutations):
    from utils import get_sorted_amino_acid_list, pairs_frequency_bar_plot
    from collections import Counter
    upload = False
    plot_key = project_name = 'outliers_two_standard_deviation_amino_acid_pairs_group_frequency' 
    x_label = 'amino acid pairs'
    y_label = 'frequency'
    top_length = 30 # indicate how many top you like to take to look at

    sorted_amino_acid_pairs_dict = {}
    for i in get_sorted_amino_acid_list():
        for j in get_sorted_amino_acid_list():
            sorted_amino_acid_pairs_dict[",".join(sorted([i[0],j[0]]))] = 0 
    
    count_mutations_by_amino_acid_pair_group = {protein:{metric:sorted_amino_acid_pairs_dict.copy() for metric in metrics} for protein in proteins}
    for protein in proteins:
        for metric, mutations in filtered_mutations[protein].items():
            for mutation in mutations:
                # update inserted animo acid count
                count_mutations_by_amino_acid_pair_group[protein][metric][",".join(sorted([mutation.aa1,mutation.aa2]))] += 1

            # sorting the pairs by frequency and pick the top many
            sorted_list = sorted(count_mutations_by_amino_acid_pair_group[protein][metric].items(),key=lambda x:x[1],reverse=True)[:top_length]     

            x = [item[0] for item in sorted_list]
            y = [item[1] for item in sorted_list]
            collect = []
            for i in x:
                collect.append(i.split(',')[0]) 
                collect.append(i.split(',')[1])
            counter = Counter(collect)
            items = sorted(counter.items(), key=lambda x:x[1],reverse=True)[:5]
            inner_data = items

            # set dummy label
            dummy_label = ""
            dummy_label += "".join([" ".join([str(x_l) + ":",str(value)+"\n"]) for x_l, value in items])
            dummy_label = dummy_label[:-1]
            
            # Create a dummy handle
            from matplotlib.lines import Line2D
            dummy_handle = Line2D([], [], color='none', label=dummy_label)
            
            plot_title = "_".join([protein,metric,plot_key])

            pairs_frequency_bar_plot(
              x,y,x_label,y_label,
              dummy_handle,
              plot_title+f"_top_{top_length}",
              upload,
              project_name,
              plot_key,
              inner_data
              )
    
def task_group_mutations_by_position_pairs(proteins,metrics,filtered_mutations):
    from utils import get_sorted_amino_acid_list, pairs_frequency_bar_plot, get_protein_chain_length

    from collections import Counter
    upload = False
    plot_key = project_name = 'outliers_two_standard_deviation_position_pairs_group_frequency' 
    x_label = 'position pairs'
    y_label = 'frequency'
    top_length = 30 # indicate how many top you like to take to look at

    


    count_mutations_by_position_pair_group = {protein:{metric:{} for metric in metrics} for protein in proteins}
    for protein in proteins:
        protein_chain_length = get_protein_chain_length(protein)
        sorted_position_pairs_dict = {",".join(sorted([str(i),str(j)])):0 for i in range(1,protein_chain_length+1+2) for j in range(1,protein_chain_length +1 + 2)}

        for metric, mutations in filtered_mutations[protein].items():
            count_mutations_by_position_pair_group[protein][metric] = sorted_position_pairs_dict.copy()
            for mutation in mutations:
                # update inserted animo acid count
                count_mutations_by_position_pair_group[protein][metric][",".join(sorted([mutation.indel1,mutation.indel2]))] += 1

            # sorting the pairs by frequency and pick the top many
            sorted_list = sorted(count_mutations_by_position_pair_group[protein][metric].items(),key=lambda x:x[1],reverse=True)[:top_length]     

            x = [item[0] for item in sorted_list]
            y = [item[1] for item in sorted_list]
            collect = []
            for i in x:
                collect.append(i.split(',')[0]) 
                collect.append(i.split(',')[1])
            counter = Counter(collect)
            items = sorted(counter.items(), key=lambda x:x[1],reverse=True)[:5]
            inner_data = items

            # set dummy label
            dummy_label = ""
            dummy_label += "".join([" ".join([str(x_l) + ":",str(value)+"\n"]) for x_l, value in items])
            dummy_label = dummy_label[:-1]
            
            # Create a dummy handle
            from matplotlib.lines import Line2D
            dummy_handle = Line2D([], [], color='none', label=dummy_label)
            
            plot_title = "_".join([protein,metric,plot_key])

            pairs_frequency_bar_plot(
                x,y,x_label,y_label,
                dummy_handle,
                plot_title+f"_top_{top_length}",
                upload,
                project_name,
                plot_key,
                inner_data
                )

class Table:
    def __init__(self, protein_name):
        self.protein_name = protein_name

        # metric in cols
        self.data = {
            'rigidity_order_parameter': [0]*10,
            'cluster_configuration_entropy': [0]*10,
            'hbond_count': [0]*10
        }
        
        # groups in rows
        self.rows = ['VS', 'S', 'M', 'L', 'VL','VSVS','SS','MM','LL','VLVL']
        
    def set_frequency(self, property_name, values):
        self.data[property_name] = values
    
    def create_template_dataframe(self):
        self.df = pd.DataFrame(self.data, index=self.rows)
        return self.df
    
    def increment_value(self, property_name, row, increment_value):
        # if property_name in self.data and row in self.rows:
        self.df.at[row, property_name] += increment_value


def task_group_mutations_by_amino_acid_size(proteins,filtered_mutations,cropped_mutations_groups,metrics):
    from group_amino_acid import get_amino_acid_size_labeled_dict
    labeled_amino_acid_size_dict = get_amino_acid_size_labeled_dict().copy()

    for protein in proteins:

        # for exhaustive mutations
        protein_template_full = Table(protein)
        df_full = protein_template_full.create_template_dataframe()
        for mutation in cropped_mutations_groups[protein]:
            aa1 = mutation.aa1
            aa2 = mutation.aa2
            size_group1 = labeled_amino_acid_size_dict[aa1]
            size_group2 = labeled_amino_acid_size_dict[aa2]
            # if size_group1 == size_group2 and size_group1 != 'S' and size_group2 != 'L':
            for metric in metrics:
                if size_group1 == size_group2:
                    protein_template_full.increment_value(metric,size_group1+size_group2,1)
                    protein_template_full.increment_value(metric,size_group1+size_group2,1)
                protein_template_full.increment_value(metric,size_group1,1)
                protein_template_full.increment_value(metric,size_group2,1)
        print(f'finish loding exhaustive mutations: {protein}')

        # for outliers
        protein_template = Table(protein)
        df = protein_template.create_template_dataframe()
        for metric, mutations in filtered_mutations[protein].items():
            for mutation in mutations:
                aa1 = mutation.aa1
                aa2 = mutation.aa2
                size_group1 = labeled_amino_acid_size_dict[aa1]
                size_group2 = labeled_amino_acid_size_dict[aa2]
                # if size_group1 == size_group2 and size_group1 != 'S' and size_group2 != 'L':
                if size_group1 == size_group2:
                    protein_template.increment_value(metric,size_group1+size_group2,1)
                    protein_template.increment_value(metric,size_group1+size_group2,1)
                protein_template.increment_value(metric,size_group1,1)
                protein_template.increment_value(metric,size_group2,1)
        
        # normalization
        protein_template_final = Table(protein)
        df_final = protein_template_final.create_template_dataframe()
        df_final = df/df_full

        # print protein name and the result
        print(protein)
        print(df_final)
        
        # clear mem
        del protein_template, df
        del protein_template_full, df_full
        del protein_template_final, df_final
    


def outlier_two_stand_deviation(cropped_mutations_groups):
    from protein import _mutants_outside_two_stdev_for_metric
    proteins = [item [0] for item in cropped_mutations_groups.items()]

    # metric list
    metrics = sorted([
        # 'size_of_largest_cluster',
        'rigidity_order_parameter', # ROP
        'cluster_configuration_entropy', # CCE
        'hbond_count' #HBond
    ])
    
    # filter the mutations to get the outliers (using 2 standard deviations) base on each metric in use
    filtered_mutations = {protein:{} for protein in proteins}
    for protein in proteins:
        for metric in metrics:
            filtered_mutations[protein][metric] = _mutants_outside_two_stdev_for_metric(cropped_mutations_groups[protein],metric)


    # task group mutations by amino acid
    # task_group_mutations_by_amino_acid(proteins,metrics,filtered_mutations)

    # task group mutations by position
    # task_group_mutations_by_position(proteins,metrics,filtered_mutations)
    
    # task group mutations by amino acid pair
    # task_group_mutations_by_amino_acid_pairs(proteins,metrics,filtered_mutations)

    # task group mutations by position pair
    # task_group_mutations_by_position_pairs(proteins,metrics,filtered_mutations)

    # task group mutations by amino acid size
    task_group_mutations_by_amino_acid_size(proteins,filtered_mutations,cropped_mutations_groups,metrics)
