#!/usr/bin/python
import glob
import sys
import math
from optparse import OptionParser
import os
import pybedtools
from pybedtools import BedTool
from os.path import isfile, join
import gzip

parser = OptionParser()

parser.add_option("--bd",
                  dest="b_input_dir",
                  help="Directory of 'b' files using wildcard '/*' notation")

parser.add_option("--bf",
                  dest="b_input_file",
                  help="Single 'b' file")

parser.add_option("--af",
				  dest="a_input_file",
				  help="Single 'a' file")

parser.add_option("--od",
                  dest="output_dir",
                  help="Output directory")

parser.add_option("--of",
                  dest="output_file",
                  help="Output file")

(options, args) = parser.parse_args()
a = pybedtools.BedTool(options.a_input_file)

if (options.output_file):
	cat = open(options.output_file, 'w')

for f in os.listdir(options.b_input_dir):
	fin = join(options.b_input_dir,f)
	fout = join(options.output_dir,f)
	b = pybedtools.BedTool(fin)
	results = a.intersect(b, wb = True).saveas('tmp')
	fread = open('tmp', 'r')
	with gzip.open(fout, 'w') as fwrite:
	    for line in fread:
			A = line.rstrip().split()
			fwrite.write('\t'.join(A[5:]) + '\n') 
			if cat:
				cat.write('\t'.join(A[5:]) + '\n')
