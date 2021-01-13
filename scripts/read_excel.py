import sys
import re
import os
import numpy as np
import pandas as pd
from optparse import OptionParser
import glob
import math

parser = OptionParser()

if len(sys.argv) < 4:
    sys.stderr.write('usage:\t' + \
                      sys.argv[0] + \
                     ' <excel ref file>' + \
                     ' <excel map file>' + \
                     ' <out file>' )
    sys.exit(1)

ref_file=sys.argv[1]
map_file=sys.argv[2]
out_file=sys.argv[3]

parser.add_option("-r",
				  "--ref_sheet",
				  dest="ref_sheet_id",
				  default=0,
				  help="Excel ref sheet name")

parser.add_option("-m",
                  "--map_sheet",
                  dest="map_sheet_id",
                  default=0,
                  help="Excel map sheet name")

parser.add_option("-c",
				  "--column",
				  dest="column_id",
				  type="int",
				  default=0,
				  help="Gene column")
	
(options, args) = parser.parse_args()

ref_data = pd.read_excel(ref_file, sheet_name = options.ref_sheet_id)
map_data = pd.read_excel(map_file, sheet_name = options.map_sheet_id)
#print ref_data
#print map_data

ref_genes = ref_data['Unnamed: 12'].tolist()
#print ref_genes
map_genes = map_data['Gene Name'].tolist()
#print map_genes
#print (ref_data.columns)
ref_coord_col = '*As the FANTOM5 consortium performed Y chromosome alignment for male and unlabeled, but not female samples, the analysis on the Y chromosome was left out.'
ref_gene_col = 'Unnamed: 12'
ref_xmale_col = 'Unnamed: 3'
ref_xfemale_col = 'Unnamed: 4'

c = 0.001

f = open(out_file, 'w')
for index, row in ref_data.iterrows():
	for find in map_genes:
		if row[ref_gene_col] == find:
			male_exp = row[ref_xmale_col] + c
			female_exp = row[ref_xfemale_col] + c
			log2fc = math.log((male_exp/female_exp), 2)
			f.write(row[ref_coord_col] + '\t' + row[ref_gene_col] + '\t' + str(log2fc) + '\n')
f.close()
