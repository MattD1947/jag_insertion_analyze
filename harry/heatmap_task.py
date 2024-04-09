import os
from os import makedirs
from pandas import read_csv as read_csv
from matplotlib.patches import Rectangle

def main():
  # feed output to feeding
  root_jag = os.path.join(os.getcwd())
  # feeding_path = os.path.join('root_jag','output')
  amino_acid_feeding_path = os.path.join(root_jag, '../result/outlier_amino_acid_pairs_frequency')
  amino_acid_feeding_list = os.listdir(amino_acid_feeding_path)
  # print(feeding_list)

  # feeding filenames
  feedings = [feeding for feeding in sorted(amino_acid_feeding_list)
  if feeding.endswith('a1_a2_frequency.txt')]

  frequency_amino_acid_data_sets = []
  for feeding in feedings:
    try:
      file_path = os.path.join(amino_acid_feeding_path, feeding)
      # print(file_path)
      frequency_amino_acid_data_sets.append(read_csv(file_path, delimiter = '\t'))
    except Exception as e:
      print(f"Error reading {feeding}: {e}")

  amino_acid_heatmap = os.path.join(root_jag,'amino_acid_heatmap')
  project_name = "rmsd_amino_acid_heatmap_by_LAOJIEWOXIANG"
  print(frequency_amino_acid_data_sets)
  # exit(1)
  [create_heatmap(data,'pair',f"{feeding[:-4].replace('_',' ')}",
                  amino_acid_heatmap,'Amino Acid','Amino Acid', False, project_name)
  for data, feeding in zip(frequency_amino_acid_data_sets, feedings)]

def get_sorted_amino_acid_list():
  sorted_amino_acid_list =sorted(
          [
          ('A',88.6 ,'Hydrophobic'),
          ('V',140.0,'Hydrophobic'),
          ('F',189.9,'Hydrophobic'),
          ('P',112.7,'Hydrophobic'),
          ('L',166.7,'Hydrophobic'),
          ('I',166.7,'Hydrophobic'),
          ('R',173.4,'Hydrophilic'),
          ('D',111.1,'Hydrophilic'),
          ('E',138.4,'Hydrophilic'),
          ('S',89   ,'Hydrophilic'),
          ('C',108.5,'Hydrophilic'),
          ('N',114.1,'Hydrophilic'),
          ('Q',143.8,'Hydrophilic'),
          ('H',153.2,'Hydrophilic'),
          ('T',116.1,'Hydrophilic,Amphipathic'),
          ('K',168.6,'Hydrophilic,Amphipathic'),
          ('Y',193.6,'Hydrophilic,Amphipathic'),
          ('M',162.9,'Amphipathic'),
          ('W',227.8,'Amphipathic'),
          ('G',60.1 ,'None')
          ],
          key=lambda x: x[1], reverse=False
      )
  return sorted_amino_acid_list
sorted_list = get_sorted_amino_acid_list()

def create_heatmap(data,pair,plot_title,output_path,x_axis_label,y_axis_label,upload,project_name):
  import matplotlib.pyplot as plt
  import seaborn as sns
  import pandas as pd
  import numpy as np

  plt.rcParams.update({'font.size': 12 * 1.8})
  sorted_amino_acid_size = [tup[0] for tup in get_sorted_amino_acid_list()]
  heatmap_data = data.pivot(index=pair, values="frequency")
  heatmap_data = heatmap_data.reindex(index=sorted_amino_acid_size[::-1], columns=sorted_amino_acid_size)
  original_index = heatmap_data.index.tolist()
  original_columns = heatmap_data.columns.tolist()
  # Combine a1, a2 and a2, a1 pairs frequency
  combined_heatmap_data = heatmap_data + heatmap_data.T

  # keep the original indexes
  combined_heatmap_data = combined_heatmap_data.reindex(index=original_index, columns=original_columns)
  print(combined_heatmap_data)

  # Create a mask for the upper right triangle
  mask_upper_right = np.triu(np.ones_like(heatmap_data, dtype=bool), k=1)
  mask_upper_left = np.fliplr(mask_upper_right)

  # Invert the mask to get the lower right triangle
  mask_lower_right = ~mask_upper_left

  # Apply the inverted mask to the combined heatmap data to get the lower right triangle
  final_heatmap_data = combined_heatmap_data.where(mask_lower_right)

  # Creating the heatmap
  plt.figure(figsize=(12, 10))  # Adjust the size as needed
  ax = sns.heatmap(final_heatmap_data, annot=False, cmap='viridis')  # 'annot=True' to display the frequency values

  # Customizations
  # plt.title(plot_title+" heatmap")
  plt.xlabel(x_axis_label)
  plt.ylabel(y_axis_label)
  # Define the rectangle starting from (0, 0), with the width and height equal to the number of columns and rows respectively
  rect = Rectangle((0, 0), width=heatmap_data.shape[1], height=heatmap_data.shape[0], fill=False, edgecolor='black', lw=2)
  ax.add_patch(rect)
  # Make color bar bolder
  cbar = ax.collections[0].colorbar
  cbar.outline.set_linewidth(1)
  plt.tight_layout()

  # Show the heatmap
  os.makedirs(output_path, exist_ok=True)  
  plt.savefig(os.path.join(output_path,plot_title.replace(' ','_')+'_heatmap.png'))

# sorted list of aa pairs based on size
def get_sorted_two_AA_lines():
    amino_acid_list = get_sorted_amino_acid_list()
    result = []
    for a1 in amino_acid_list:
        for a2 in amino_acid_list:
            result.append(((a1[0],a2[0]),a1[1]+a2[1],a1[1],a1[2],a2[1],a2[2]))
    result = sorted(result,key=lambda x: x[1], reverse=True)
    return result

if __name__ == "__main__":
  main()

