def get_amino_acid_size_labeled_dict():
   VS = {'VS':'AGS'}
   S = {'S':'NDCPT'}
   M = {'M':'QEHV'}
   L = {'L':'RILKM'}
   VL = {'VL':'FWY'}
   amino_acid_size_dict = {}
   collect = [
      VS,
      S,
      M,
      L,
      VL
   ]

   for i in collect:
      for group_name,item in i.items():
         for j in item:
               amino_acid_size_dict[j] = group_name
   return amino_acid_size_dict
