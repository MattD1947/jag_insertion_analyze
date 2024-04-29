"""
How Pairs of Insertion Mutations Impact Protein Structure: An Exhaustive Computational Study

2024

Changrui Li, Yang Zheng, and Filip Jagodzinski
"""
class Protein:
    def __init__(self, name):
        # protein name
        self.name = name
        
        # mutations list
        self.mutations = []
        
    def __repr__(self):
        return f"<Protein name={self.name}, number of mutations={len(self.mutations)}>"

class Mutation:
    def __init__(self, indel1, aa1, indel2, aa2, hbond_count, size_of_largest_cluster, rigidity_order_parameter, cluster_configuration_entropy):
        
        # the first position
        self.indel1 = indel1
        
        # the first inserted amino acid
        self.aa1 = aa1
        
        # the second position
        self.indel2 = indel2
        
        # the second inserted amino acid
        self.aa2 = aa2
        
        # metrics
        self.hbond_count = hbond_count
        self.size_of_largest_cluster = size_of_largest_cluster
        self.rigidity_order_parameter = rigidity_order_parameter
        self.cluster_configuration_entropy = cluster_configuration_entropy

    def __repr__(self):
        return f"<Mutation indel1={self.indel1}, aa1={self.aa1}, indel2={self.indel2}, aa2={self.aa2}, hbond={self.hbond_count}, size of largest cluster={self.size_of_largest_cluster}, rigidity order parameter={self.rigidity_order_parameter}, cluster configuration entropy={self.cluster_configuration_entropy}>"