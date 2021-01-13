import sys
import re
import os
import os.path
import requests
import wget

if len(sys.argv) != 3:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <ENCODE download file>' + \
                     ' <out dir>' )
    sys.exit(1)

download_file=sys.argv[1]
out_dir=sys.argv[2]
encode_path = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/'

""" Define a function that downloads files """
def downloadfiles(infile,path,outdir):
	f1 = open(infile, 'r') 
	for line in f1:
		A = line.rstrip().split()
		url = path + A[1]
		fname = outdir + '/' + A[1].split('/')[1].split('.')[0] + '.bed.gz'
		wget.download(url,fname)
		
downloadfiles(download_file, encode_path, out_dir)

