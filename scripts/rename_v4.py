import sys
import re
import os

if len(sys.argv) != 6:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <name2library file>' + \
                     ' <expression count matrix file>' + \
                     ' <out dir>' + \
                     ' <liftOver dir>' + \
                     ' <escapee file>\n')
    sys.exit(1)

name2library_file=sys.argv[1]
expression_count_matrix_file=sys.argv[2]
out_dir=sys.argv[3]
liftOver_path = sys.argv[4]
escapee_file = sys.argv[5]

chain = "hg19ToHg38.over.chain.gz"
tmp0 = "tmp0.bed"
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

""" Define a function that filters escapee TSSs """
def filterExpression(exp,escapee_file):
	f = open(tmp0, 'a+')
	regex_to_gene = {}
		
	for l in open(escapee_file, 'r'):
		A = l.rstrip().split()
		B = re.split(r'_', A[0])
		if len(B) == 4:
			str = B[3]
		else:
			str = "_".join(B)
		regex_to_gene[str] = A[1]

	for regex in regex_to_gene.keys():
		regexp = re.compile(regex)
		with open(exp) as fin:
			for line in fin:
				if regexp.search(line):
					f.write(regex_to_gene[regex] + '\t' + line.rstrip() + '\n')

	f.close()

""" Define a function that makes header array and splits up file into 2 temp files according to genome assembly """
def parseExpression(exp,exptmp,header,tmp1,tmp2):
	line = 0
	t1 = open(tmp1,'w')
	t2 = open(tmp2,'w')
	tmpexpin = open(exptmp, 'r')
	with open(exp) as f1:
		data = f1.readlines()[1837]
		A = data.rstrip().split()
		header = A[7:]
		for i, x in enumerate(header):
			match = re.search(r'\CNh\w+', x)
			header[i] = match.group()
	with open(exptmp) as f2:
		for l in f2:
			A = l.rstrip().split()
			print A[1]
			id = A[0]
			#B = re.split('\.{2}|[:;,]+', A[1])
			#print len(B)
			#print B[0]
			#print B[1]
			#print B[3]
			#print B[4]
			#if re.match(r'hg19', l):
			#	print >> t1, chr + '\t' + start + '\t' + end + '\t' + id + '\t' + strand + '\t' + '\t'.join(A[-1829:])
			#else:
			#	print >> t2, chr + '\t' + start + '\t' + end + '\t' + id + '\t' + strand + '\t' + '\t'.join(A[-1829:])

	t1.close()
	t2.close()
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
    

#filterExpression(expression_count_matrix_file, escapee_file)
header = parseExpression(expression_count_matrix_file, tmp0, header, hg19tmp, hg38tmp)
#splitBed(hg19tmp, tmp1, tmp2)
#doliftOver(liftOver_path, chain, tmp1)
#liftOver_file = tmp1 + ".success"
#mergeliftOver(liftOver_file, tmp2, hg38tmp)
#filesExpression(hg38tmp, names, header, files, out_dir)

#for f in os.listdir(os.curdir):
#	if re.search(r'tmp', f):
#		os.remove(f)
