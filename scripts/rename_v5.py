import sys
import re
import os
import numpy as np
from optparse import OptionParser
import glob

parser = OptionParser()

if len(sys.argv) < 3:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <escapee file>' + \
                     ' <out dir>' )
    sys.exit(1)

escapee_file=sys.argv[1]
out_dir=sys.argv[2]

parser.add_option("-a",
                  "--all",
                  dest="all_file",
                  help="Concatenated gene file name")

(options, args) = parser.parse_args()
tmp0 = "tmp0.bed"

""" Define a function that makes header array and splits up file into 2 temp files according to genome assembly """
def parseExpression(esc,header,tmp):
	temp = open(tmp0,'w')
	f1 = open(esc, 'r') 
	for line in f1:
		A = line.rstrip().split()
		header.append(A[1]) if A[1] not in header else header
		[chr, start, stop, strand] = re.split('\.{2}|[:,]+', A[0])
		start = str(int(start) - 500)
		stop = str(int(stop) + 500)
		print >> temp, chr + '\t' + start + '\t' + stop + '\t' + A[1] + '\t' + strand + '\t' + A[2]

	temp.close()
	return header

""" Define a function that creates one bed file per sampe in expression header """
def filesExpression(infile,header,files,dir):
	file = open(infile, 'r')
	for l in file:
		A = l.rstrip().split()
		fname = dir + '/' + A[3] + '.bed'
		if fname not in files:
			files.append(fname)
			with open(fname, 'w') as fwrite:
				fwrite.write(l)

""" Define a function that filters gene file for TSS with largest differential expression """
def filterExpression(dir):
	for f in os.listdir(dir):
		max = np.inf
		lwrite = ''
		fname = dir + '/' + f
		fread = open(fname, 'r')
		for l in fread:
			A = l.rstrip().split()
			if float(A[5]) < float(max):
				lwrite = '\t'.join(A[:-1])
				max = A[5]
		fwrite = open(fname, 'w')
		fwrite.write(lwrite)

""" Define a function that cats gene files """
def joinGenes(dir,outfile):
	fwrite = open(outfile, 'w')
	for f in os.listdir(dir):
		fname = dir + '/' + f
		fread = open(fname, 'r')
		for l in fread:
			fwrite.write(l + '\n')

files = []
header = []

header = parseExpression(escapee_file, header, tmp0)
filesExpression(tmp0, header, files, out_dir)

filterExpression(out_dir)
if(options.all_file):
	joinGenes(out_dir,options.all_file)

for f in os.listdir(os.curdir):
	if re.search(r'tmp', f):
		os.remove(f)
