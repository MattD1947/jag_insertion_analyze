sorted_amino_acid_list =sorted(
    [
    ('A',88.6 ,'Hydrophobic'),
    ('V',99.13,'Hydrophobic'),
    ('F',189.9,'Hydrophobic'),
    ('P',112.7,'Hydrophobic'),
    ('L',166.7,'Hydrophobic'),
    ('I',166.7,'Hydrophobic'),
    ('R',173.4,'Hydrophilic'),
    ('D',111.1,'Hydrophilic'),
    ('E',138.4,'Hydrophilic'),
    ('S',89.0 ,'Hydrophilic'),  
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
    key=lambda x: x[1], reverse=False #  sorted by size in ascending order
)


# smallest and largest amino acid in volume or say size
smallest_size = sorted_amino_acid_list[0][1]
largest_size = sorted_amino_acid_list[-1][1]

