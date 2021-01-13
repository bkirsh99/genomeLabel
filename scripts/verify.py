#!/usr/bin/python
import glob
import sys
from optparse import OptionParser
import os
import csv
import numpy as np
import re

parser = OptionParser()

parser.add_option("-i",
                  "--index_dir",
                  dest="index_dir",
                  help="Index directory")

parser.add_option("-q",
                  "--query_dir",
                  dest="query_dir",
                  help="Query directory")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-m",
				  "--map",
				  dest="map_file",
				  help="ENCODE url to factor file name")

(options, args) = parser.parse_args()
map = {}
results  = {}
arr = []
arr.append(['Antibody', 'Expected', 'Observed'])

for file_name in glob.glob(options.index_dir):
	name = file_name.split('/')[1].split('.')[0]
	arr[0].append(name)

with open(options.map_file, 'r') as f0:
	reader = csv.reader(f0)
	for row in reader:
		A = row[0].rstrip().split()
		factor_name = A[0]
		url = A[1]
		expected_overlap = A[2]
		arr.append([url, expected_overlap])
		map[url] = factor_name

for list in arr[1:]:
	for match in map.keys():
		for file_name in glob.glob(options.query_dir):
			regex = file_name.split('/')[-1].split('.')[0]
			if re.search(regex, match) and list[0] == match:
				list[0] = map[match]
				f = open(file_name, 'r')
				next(f)
				for l in f:
					A = l.strip().split('\t')
					list.append(A[2])

for list in arr[1:]:
	total = 0
	for observed_overlap in range(3, len(list)): 
		if int(list[observed_overlap]) > 0:
			total = total + 1
	list.insert(2, total)

with open(options.output_file, 'w') as file:
    writer = csv.writer(file, delimiter='\t', lineterminator='\n')
    writer.writerows(arr)
