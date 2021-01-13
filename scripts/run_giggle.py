import sys
import re
import os
import numpy as np
from optparse import OptionParser
import glob
import subprocess

parser = OptionParser()

if len(sys.argv) < 2:
	sys.stderr.write('usage:\t' + \
                      sys.argv[0] + \
                     ' <input directory>' + \
                     ' <output directory>')
	sys.exit(1)

rootDir=sys.argv[1]
resultsDir=sys.argv[2]

stats = ('odds','sig','combo')

parser.add_option("-i",
                  "--index",
                  dest="index_label",
                  help="Annotation of index file(s)")

parser.add_option("-q",
				  "--query",
				  dest="query_label",
				  help="Annotation of query file(s)")

(options, args) = parser.parse_args()
indexLabel = options.index_label
queryLabel = options.query_label

""" Define a function to traverse directory tree """
def traverseRoot(rootDir):
	global indexDirs
	global queryDirs
	for dirName, subdirList, fileList in os.walk(rootDir,topdown=False):
		if re.search('/'+indexLabel+'$',dirName):
			indexDirs = doIndex(dirName)
			#print indexDirs
			#traverseSub(indexName,fileList)
		elif re.search('/'+queryLabel+'$',dirName):
			queryDirs = doIndex(dirName)
			#print queryDirs
			#traverseSub(queryName,fileList)

""" Define a function to traverse index and query directory tree """
def traverseSub(subDir,fList):
	for fname in fList:
		file = os.path.join(subDir,fname)
		print file

""" Define a function to concatenate files in directory """
def concatFiles(inDir,outFile):
	print inDir
	fwrite = open(outFile, 'w')
	for f in os.listdir(inDir):
		#fname = dir + '/' + f
		print f
		#with open(f, 'r') as fread:
			#fwrite.write(fread.read() + '\n')

""" Define a function to run GIGGLE index """
def doIndex(inDir):
	splitDir = inDir.rstrip().split('/',2)[-1].replace('/','_') + '_splitsort'
	splitIndexDir = splitDir + '_b'
	if not os.path.exists(splitDir):
		try:
			os.mkdir(splitDir)
			sort_command = 'bash $GIGGLE_ROOT/scripts/sort_bed "%s/*.bed" %s/ > /dev/null' % (inDir,splitDir)
			os.system(sort_command)
			index_command = '$GIGGLE_ROOT/bin/giggle index -i "%s/*gz" -o %s -f -s 2> /dev/null' % (splitDir,splitIndexDir)
			os.system(index_command)
		except OSError:
			print ("Creation of the directory %s failed" % splitDir)

	return splitDir, splitIndexDir
	
""" Define a function to run GIGGLE search """
def doSearch(indexDir,queryDir,outDir):
	if not os.path.exists(outDir):
		try:
			os.mkdir(outDir)
			search_command = "ls %s | ./gargs -p 10 '$GIGGLE_ROOT/bin/giggle search -i %s -q %s/{} -s > %s/{}.results'" % (queryDir,indexDir,queryDir,outDir)
			os.system(search_command)
		except OSError:
			print ("Creation of the directory %s failed" % outDir)

""" Define a function to make heat map of results """
def plotHeatmap(inDir):
	for s in stats:
		figName = inDir + '_' + s + '.pdf'
		plot_command = "./plot_general.py -i '%s/*' -o %s --stat %s" % (inDir,figName,s)
		os.system(plot_command)

traverseRoot(rootDir)
i = indexDirs[1]
q = queryDirs[0]
doSearch(i,q,resultsDir)
#plotHeatmap(resultsDir) 
#concatFiles(indexName,'try.bed')
