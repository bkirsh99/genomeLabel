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
		factor = A[0]
		#factor = A[0].replace('(','').replace(')','').replace('_','-')
		url = path + A[1]
		fname = outdir + '/' + factor + '_sort_peaks.narrowPeak.bed.gz'
		if not os.path.exists(fname):
			print fname
			#wget.download(url,fname)
			#r = requests.get(url, allow_redirects=True)
			#open(fname, 'wb').write(r.content)
		else:
			print fname

downloadfiles(download_file, encode_path, out_dir)

