# How Pairs of Insertion Mutations Impact Protein Structure: An Exhaustive Computational Study
Changrui Li, Yang Zheng, Filip Jagodzinski | Western Washington University, Department of Computer Science, 2024

## Dependency
````
pip install requirements.txt
````

## Data
Extract `data/saved_classes_data.rar` file locally under `data` folder so the script can run.

For linux, you can use the command below if you have `unrar`.
````
cd data
unrar x saved_classes_data.rar ./
cd ..
````

## Configuration
In `config/main.yaml` file, there are string of proteins and metrics the study used.

## Data Structure
Find at `src/protein.py`.

## Sample
Get into source folder.
```bash
cd src
```
If you are looking for a summary among the three proteins, run the commands below.

```bash
python3 main.py
```

If you like to look at a specific mutation of the a protein as example, for protein 1hhp, insert the first position at 1 with amino acid Proline (P), and insert the second position at 22 with alanine (A), you can run the command below.

```bash
python3 main.py 1hhp 1 P 22 A
```

And you would expect getting result like below.
<pre style="white-space: pre-wrap: break-word;">
<code>
[<Mutation indel1=1, aa1=P, indel2=22, aa2=A, hbond=10, size of largest cluster=126.0, rigidity order parameter=0.0616137, cluster configuration entropy=2.64985>]
</code>
</pre>

## Author Contributions Statement
F.J. conceived the initial experiments. C.L. and Y.Z. conducted the experiments and analyzed the results. C.L.,Y.Z., and F.J.
wrote the manuscript.

## Acknowledgements
The authors thank Zach McGrew in the CS Support Group for his expertise in maintaining the compute cluster. This work was supported in part by funds from the National Science Foundation (NSF: #2031283).

# Proteins and Metrics

proteins: 
1crn, 1csp, 1hhp

metrics: 
cluster_configuration_entropy,
hbond_count,
size_of_largest_cluster,
rigidity_order_parameter
