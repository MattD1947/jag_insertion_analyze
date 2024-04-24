"""
Identifying And Comparison Analyze Impactful Pairs of Insertion Mutations in Proteins

2024

Changrui Li and Yang Zheng
"""

# Standard libary imports
import os

# Third party imports
import yaml, json, pickle, statistics
# wandb
import matplotlib.pyplot as plt
import numpy as np

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

# create sorted amino acid lsit
def get_sorted_amino_acid_list():
  sorted_amino_acid_list =sorted(
          [
          ('A',88.6 ,'Hydrophobic'),
          ('V',140.00,'Hydrophobic'),
          ('F',189.9,'Hydrophobic'),
          ('P',112.7,'Hydrophobic'),
          ('L',166.7,'Hydrophobic'),
          ('I',166.7,'Hydrophobic'),
          ('R',173.4,'Hydrophilic'),
          ('D',111.1,'Hydrophilic'),
          ('E',138.4,'Hydrophilic'),
          ('S',89.0,'Hydrophilic'),
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

# read protein chain A length
def get_protein_chain_length(protein):
    link = "https://files.rcsb.org/header/"+protein+".pdb"
    import requests
    lines = requests.get(link).text
    lines = lines.split('\n')
    length = 0
    for line in lines:
        if line[:5] == 'DBREF':
            length = int(line.split()[4])
            break
    return length

# return positions that indels on alpha helix secondary structrue
def get_on_alpha_helix_positions(protein_name):
  helix_1hhp = [ i - 1 for i in range(86,94+1)]
  helix_1crn = [i - 1 for i in range(7,19+1)] + [j - 1 for j in range(23,30+1)]
  helix_1csp = [30 - 1,31 - 1,32 - 1]
  data = {
      '1hhp':helix_1hhp,
      '1crn':helix_1crn,
      '1csp':helix_1csp
  }
  return data[protein_name]

# upload to wandb
# def upload_to_wandb(plt,project_name,plot_key,image_title):
#     wandb.init(
#         entity='jag_research',
#         project=project_name,
#         name=image_title
#         )
#     wandb.log({plot_key:wandb.Image(plt)})
#     wandb.finish()

def single_frequency_bar_plot(x,y,x_label,y_label,dummy_handle,plot_title,upload,project_name,plot_key):
  # fontsize
  fontsize =14*1.8
  x_font_size = min(fontsize,(1 + (3.87 - np.log(len(x)))/2)*fontsize)
  # print(x_font_size)

  multiple = 1.0


  plt.figure(figsize=(12*multiple,8*multiple))
  plt.bar(x,y,color='green',edgecolor='black')

  if 'position' in plot_title:
    protein_name = plot_title[:4]
    helix_position = get_on_alpha_helix_positions(protein_name)
    plt.bar(helix_position,[max(y)]*len(helix_position),color='red',edgecolor='red',alpha=0.3)

  # plt.title(plot_title)
  plt.xlabel(x_label,fontsize=fontsize)
  plt.ylabel(y_label,fontsize=fontsize)
  # inner = max(y) - min(y)
  # bottom = min(y) - int(inner/10)
  # top =  max(y) + int(inner/10)
  # plt.ylim(bottom,top)
  # print(plot_title,bottom,top)

  plt.xticks( fontsize = fontsize) if x[0].isalpha() else  plt.xticks(x[::10],fontsize=fontsize)
  plt.yticks(fontsize = fontsize)
  plt.grid(False)

  # Append the dummy handle and the new label
  handles, labels = plt.gca().get_legend_handles_labels()
  handles.append(dummy_handle)
  labels.append(dummy_handle.get_label())


  # Set the updated legend
  plt.legend(handles=handles, labels=labels,fontsize=fontsize, handlelength=0, handletextpad=0)
  plt.tight_layout()
  # plt.show()
  # plt.close()

  if not upload:
    # plt.show()
    folder_path = os.path.join('../plots_results',project_name)
    from os import makedirs
    makedirs(folder_path,exist_ok=True)
    plt.savefig(os.path.join(folder_path,plot_title+'.png'),dpi=400)
  else:
    # upload_to_wandb(plt,project_name,plot_key,plot_title)
    pass
  plt.close()

def pairs_frequency_bar_plot(x,y,x_label,y_label,dummy_handle,plot_title,upload,project_name,plot_key,inner_data):
  # fontsize
  fontsize =14*1.8

  # plt.figure(figsize=(13,8))
  # plt.bar(x,y,color='green',edgecolor='black')
  # # plt.title(plot_title)
  # plt.xlabel(x_label,fontsize=fontsize)
  # plt.ylabel(y_label,fontsize=fontsize)
  # plt.xticks(rotation=90)
  # plt.grid(False)

  # Append the dummy handle and the new label

  from mpl_toolkits.axes_grid1.inset_locator import inset_axes
  fig,ax = plt.subplots(figsize=(13,8))
  bars_main = ax.bar(x,y,color='green',edgecolor='black')
  ax.set_xlabel(x_label,fontsize=fontsize)
  ax.set_ylabel(y_label,fontsize=fontsize)
  ax.set_xticklabels(x,rotation=90,fontsize=fontsize)
  bottom = min(y) - 5
  top = max(y) + 3
  ax.set_ylim(bottom,top)
  y_stick_labels = ax.get_yticklabels()
  ax.set_yticklabels(y_stick_labels, rotation=0,fontsize=fontsize)
  ax.grid(False)
  plt.tight_layout()

  # handles, labels = plt.gca().get_legend_handles_labels()
  # handles.append(dummy_handle)
  # labels.append(dummy_handle.get_label())
  # # Set the updated legend
  # plt.legend(handles=handles, labels=labels,fontsize=fontsize, handlelength=0, handletextpad=0)



  axins = inset_axes(ax,width="40%",height="55%",loc=1,borderpad=1)

  axins.patch.set_facecolor('white')
  axins.patch.set_alpha(1)

  # axins.set_visible(True)
  inner_x =[item[0] for item in inner_data]
  inner_y =[item[1] for item in inner_data]
  bars_inset = axins.bar(inner_x,inner_y,color="red")

  axins.set_ylim(min(inner_y) - 2,max(inner_y) + 2)


  for bar in bars_inset:
    yval = bar.get_height()
    axins.text(bar.get_x() + bar.get_width()/2, yval + 0.1, round(yval, 1), ha='center', va='bottom',fontsize=fontsize-6)
  # for label in axins.get_xticklabels():
  #   label.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='none'))
  axins.set_xticklabels(axins.get_xticklabels(), fontsize=fontsize-6)
  # for label in axins.get_yticklabels():
  #   label.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='none'))
  axins.set_yticklabels(axins.get_yticklabels(), fontsize=fontsize-6)


  if not upload:
    # plt.show()
    folder_path = os.path.join('../plots_results',project_name)
    from os import makedirs
    makedirs(folder_path,exist_ok=True)
    plt.savefig(os.path.join(folder_path,plot_title+'.png'),dpi=400)
  # else:
  #   upload_to_wandb(plt,project_name,plot_key,plot_title)
  #   pass
  plt.close()