import sys
import re
import os

if len(sys.argv) != 5:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <name2library file>' + \
                     ' <expression count matrix file>' + \
                     ' <out dir>' + \
                     ' <liftOver dir>\n')
    sys.exit(1)

name2library_file=sys.argv[1]
expression_count_matrix_file=sys.argv[2]
out_dir=sys.argv[3]
liftOver_path = sys.argv[4]

chain = "hg19ToHg38.over.chain.gz"
tmp1 = "tmp1.bed"
tmp2 = "tmp2.bed"
hg19tmp = "tmphg19.bed"
hg38tmp = "tmphg38.bed"

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
def parseExpression(exp,header,tmp1,tmp2):
	line = 0
	t1 = open(tmp1,'w')
	t2 = open(tmp2,'w')
	expin = open(exp,'r')
	for l in expin:
		line += 1
		A = l.rstrip().split()
		if line == 1838:
			header = A[7:]
			for i, x in enumerate(header):
				match = re.search(r'\CNh\w+', x)
				header[i] = match.group()
			0/1
		elif line >= 1840:
			[genome, chr, start, end, strand, id] = re.split(r'\.{2}|[:;,]+', A[0])
			id = A[1]
			if re.match(r'hg19', l):
				print >> t1, chr + '\t' + start + '\t' + end + '\t' + id + '\t' + strand + '\t' + '\t'.join(A[-1829:])
			else:
				print >> t2, chr + '\t' + start + '\t' + end + '\t' + id + '\t' + strand + '\t' + '\t'.join(A[-1829:])
	t1.close()
	t2.close()
	expin.close()
	return header

""" Define a function that creates one bed file per sampe in expression header """
def filesExpression(fparsed,names,header,files,dir):
	file = open(fparsed, 'r')
	for l in file:
		A = l.rstrip().split()
		k = 0
		for a in A[-1829:]:
			if a != '0':
				fname = dir + '/' + names[header[k]] + '.bed'
				if fname not in files:
					files[fname] = open(fname, 'w')
				fstr = "\t".join(A[0:4]) + '\t' + a + '\t' + A[4] + '\n'
				files[fname].write(fstr)
			k+=1

files = {}
names = {}
header = []

for l in open(name2library_file, 'r'):
    A = l.rstrip().split('\t')
    names[A[1]] = A[0]

#header = parseExpression(expression_count_matrix_file, header, hg19tmp, hg38tmp)
#splitBed(hg19tmp, tmp1, tmp2)
#doliftOver(liftOver_path, chain, tmp1)
#liftOver_file = tmp1 + ".success"
#mergeliftOver(liftOver_file, tmp2, hg38tmp)
#filesExpression(hg38tmp, names, header, files, out_dir)
#print os.getcwd()

for f in os.listdir(os.curdir):
	if re.search(r'tmp', f):
		os.remove(f)
