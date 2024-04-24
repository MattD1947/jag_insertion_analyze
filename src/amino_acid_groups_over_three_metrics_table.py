# from group_amino_acid import VS,S,M,L,VL
import pandas as pd
from pandas import read_csv

def concat_data(rmsd_data,cat1,cat2,output_root,base_name):
  # set df cat 1 and 2
  # print(rmsd_data.columns)
  df_cat1 = rmsd_data[['protein',cat1]]
  df_cat2 = rmsd_data[['protein',cat2]]

  # set new column name
  new_column = 'amino_acid' if cat1 == 'a1' and cat2 =='a2' else 'position'

  # reset colomns
  df_cat1.columns = ['protein',new_column]
  df_cat2.columns = ['protein',new_column]


  # concat df cat 1 and 2 and name it df_cat_amino_acid
  df_cat = pd.concat([df_cat1,df_cat2],ignore_index=True)

  # Group by 'amino_acid' and count the frequency
  # print(df_cat_amino_acid.keys)
  group_cat = df_cat.groupby(['protein',new_column]).size().reset_index(name='frequency') # jjj

  # Sort by frequency in descending order
  group_cat_sorted = group_cat.sort_values(by='frequency', ascending=False)
  print(group_cat_sorted)

  # Save the result to a text file
#   output_file_amino_acid_updated = os.path.join(output_root,f'{base_name}_{new_column}_frequency.txt')
#   group_cat_sorted.to_csv(output_file_amino_acid_updated, index=False, sep='\t')
#   # print()


def amino_acid_groups_over_three_metrics_table():
    outlier_rmsd_root = 'data/outlier_rmsd'
    import os
    data = {}
    for filename in os.listdir(outlier_rmsd_root):
        protein_name = filename[:4]
        if data.get(protein_name) == None:
            data[protein_name] = {}
        metric_name = filename[5:-len('_outlier_rmsd.txt')]
        
        rmsd_data_sets = read_csv(os.path.join(outlier_rmsd_root,filename))

        data[protein_name][metric_name] = rmsd_data_sets
    tmp_data = data['1hhp']['rigidity_order_parameter']
    # group_data setting
    concat_data(rmsd_data=tmp_data,cat1='a1',cat2='a2',output_root='.',base_name='base_name')


    # for file,rmsd_data in zip(filenames,rmsd_data_sets):
    #     base_name = file[:-len('_rmsd.txt')]
    #     aa_pairs_frequency(rmsd_data,output_amino_acid_root,base_name)
        # pp_pairs_frequency(rmsd_data,output_position_root,base_name)
    pass