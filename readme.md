# How Pairs of Insertion Mutations Impact Protein Structure: An Exhaustive Computational Study
Changrui Li, Yang Zheng, Filip Jagodzinski | Western Washington University, Department of Computer Science, 2024

## Dependency
```
```bash
pip install requirements.txt
```
```

## Data
Extract `data/saved_classes_data.rar` file locally under `data` folder so the script can run.


## Configuration
In `main.yaml` file, there are proteins and metrics the study used.

## Data Structure
Find at `src/protein.py`.

## Sample
Run 

```
```bash
cd src
python3 main.py
```
```

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
