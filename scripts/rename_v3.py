import sys
import re
import os

if len(sys.argv) != 4:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <escapee file>' + \
                     ' <out dir>' + \
                     ' <liftOver dir>')
    sys.exit(1)

escapee_file=sys.argv[1]
out_dir=sys.argv[2]
liftOver_path = sys.argv[3]

chain = "hg19ToHg38.over.chain.gz"
tmp0 = "tmp0.bed"

""" Define a function that implements liftOver """
def doliftOver(liftOver,chain,infile):
	cmd = " ".join([liftOver,infile,chain,infile + ".success",infile + ".failure"])
	os.system(cmd)

""" Define a function that splits up bed file into 2 temp files """
def splitBed(bed,tmp1,tmp2):
	t1 = open(tmp1,'w')
	t2 = open(tmp2,'w')
	bedin = open(bed,'r')
	for l in bedin:
		e = l.strip().split("\t")
		print >> t1, "\t".join(e[0:3])
		print >> t2, "\t".join(e[3:])
	t1.close()
	t2.close()
	bedin.close()

""" Define a function that merges liftOver """
def mergeliftOver(f1,f2,outputfile):
	with open(f1, 'r') as fh1, open(f2, 'r') as fh2, open(outputfile, 'a+') as o:
	    lines_f1 = fh1.readlines()
	    lines_f2 = fh2.readlines()
	    for i in range(len(lines_f1)):
	    	o.write(lines_f1[i].strip() + '\t' + lines_f2[i].strip())
	fh1.close()
	fh2.close()
	o.close()

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
		print >> temp, chr + '\t' + start + '\t' + stop + '\t' + A[1] + '\t' + strand

	temp.close()
	return header

""" Define a function that creates one bed file per sampe in expression header """
def filesExpression(infile,header,files,dir):
	file = open(infile, 'r')
	for l in file:
		A = l.rstrip().split()
		fname = dir + '/' + A[3] + '.bed'
		if fname not in files:
			files[fname] = open(fname, 'w')
		files[fname].write(l)

files = {}
header = []


header = parseExpression(escapee_file, header, tmp0)
#doliftOver(liftOver_path, chain, tmp0)
#liftOver_file = tmp0 + ".success"
filesExpression(tmp0, header, files, out_dir)

for f in os.listdir(os.curdir):
	if re.search(r'tmp', f):
		os.remove(f)
