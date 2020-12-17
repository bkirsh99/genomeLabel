#!/usr/bin/python
# -*- coding: future_fstrings -*-
#!pip install brewer2mpl

import glob
import sys
import math
from optparse import OptionParser
import os
import io
from scipy.cluster.hierarchy import inconsistent, fcluster,ward
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns

parser = OptionParser()

parser.add_option("-i",
                  "--input",
                  dest="input",
                  help="Input results directory or file")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("--stat",
                  dest="stat",
                  default="combo_score",
                  help="Stat to describe (query, filename, size, overlaps, odds_ratio, fishers_two_tail, fishers_left_tail, fishers_right_tail, combo_score) (Default: combo_score)")

parser.add_option("--x_stat",
				  dest="x_stat",
				  help="X-axis stat")

parser.add_option("--y_stat",
                 dest="y_stat",
				 help="Y-axis stat")

parser.add_option("--x_size",
				  dest="x_size",
                  type="int",
                  default=10,
                  help="Figure x size (Default 10)")

parser.add_option("--y_size",
                  dest="y_size",
                  type="int",
                  default=30,
                  help="Figure x size (Default 30)")

parser.add_option("--by",
				  dest="group_by",
                  default="query",
                  help="Output by query or index annotations")

parser.add_option("--no_labels",
                  dest="no_labels",
                  action="store_true",
                  default=False,
                  help="Do not label x and y axis");

parser.add_option("--no_xlabels",
                  dest="no_xlabels",
                  action="store_true",
                  default=False,
                  help="Do not label x axis");

parser.add_option("--no_ylabels",
                  dest="no_ylabels",
                  action="store_true",
                  default=False,
                  help="Do not label y axis");

parser.add_option("--row_names",
                  dest="row_map",
                  help="File to map row names")

parser.add_option("--col_names",
				  dest="col_map",
				  help="File to map column names")

parser.add_option("-n",
				  dest="n",
				  type="int",
				  default=5,
				  help="Top n to display (Default: 5)")

(options, args) = parser.parse_args()

if options.stat not in ['query', 'filename', 'size', 'overlaps', 'odds_ratio', 'fishers_two_tail', 'fishers_left_tail', 'fishers_right_tail', 'combo_score']:
    parser.error('Stat "' + options.stat + '" not supported')

def set_plot_param(arr):
	(title,font,x_fig,y_fig,label) = arr
	params = {'axes.titlesize': title,
	          'legend.fontsize': font,
	          'figure.figsize': (x_fig, y_fig),
	          'axes.labelsize': label,
	          'xtick.labelsize': label,
	          'ytick.labelsize': label,
	          'figure.titlesize': title}
	plt.rcParams.update(params)
	plt.style.use('seaborn-whitegrid')
	sns.set_style("white")

def make_scatter(df,name,a,b,params):
	categories = np.unique(df[name])
	colors = [plt.cm.tab10(i/float(len(categories)-1)) for i in range(len(categories))]
	set_plot_param(params)
	#fig = plt.figure(figsize=(16,10),dpi=80,facecolor='w',edgecolor='k')
	for i,category in enumerate(categories):
		plt.scatter(a,b,data=df.loc[df[name]==category, :],c=colors[i],label=str(category),edgecolors='black',linewidth=.5)
	plt.xlabel(a)
	plt.ylabel(b)
	plt.savefig('testlala.pdf', format='pdf')

def make_barplot(df,name,stat,params):
	df[[stat,name]].groupby(name).apply(lambda x: x.mean())
	df.sort_values(stat,inplace=True)
	df.reset_index(inplace=True)
	fig,ax=plt.subplots()
	set_plot_param(params)
	#plt.bar(df[
	for i,stat in enumerate(df[stat]):
		ax.text(i,stat+0.5,round(stat,1),horizontalalignment="center")
	plt.savefig('testlala2.pdf',format='pdf')

def parse_in(path):
	results = ''
	ith_file = 0
	files = []

	if os.path.isdir(path):
		print "is directory"
		for file in glob.glob(os.path.join(path, '*')):
			(results, ith_file, files) = filter_cluster(file,results,ith_file,files)
	elif os.path.isfile(path):
		print "is file"
		(results, ith_file, files) = filter_cluster(path,results,ith_file,files)

	res_table = pd.read_csv(io.StringIO(results),sep="\t",header=None)
	res_table.columns = ['query', 'filename', 'size', 'overlaps', 'odds_ratio', 'fishers_two_tail', 'fishers_left_tail', 'fishers_rigth_tail', 'combo_score']
	return res_table

def filter_cluster(file,str,n,arr):
	if os.stat(file).st_size != 0:
		with open(file, 'r') as f:
			result = f.read().decode('utf-8').replace('\t\n', '\n')
	n += 1
	arr += [file]
	str = str + file + '\t' + (f'\n{file}\t').join(result.split("\n")[1:-1]) + "\n"
	return (str, n, arr)
	
def scoreDistribution(df,figname):
	f = plt.figure()
	n, bins, patches = plt.hist(df.combo_score.values,color="#0504aa",bins=1600,alpha=0.7,rwidth=0.85,density=True)
	plt.xlim(-10,10)
	plt.xlabel("Combo Score", fontsize=12)
	plt.xticks(fontsize=10)
	plt.yticks(fontsize=10)
	plt.ylabel("Probability",fontsize=12)
	plt.title("Combo Score Distribution",fontsize=15)

	f.savefig(figname,dpi=None,facecolor="w",edgecolor="w",orientation="portrait",
			  papertype=None,format=None,transparent=False,bbox_inches='tight',pad_inches=0.8,
			  frameon=None,metadata=None,figsize=[6.4,4.8])

def rename_labels(df):
	df['query'] = df['query'].apply(lambda x: os.path.basename(x))
	df['filename'] = df['filename'].apply(lambda x: os.path.basename(x))
	df = df.replace(regex=r'.bed.gz.*$',value='')
	return df

def group_dataframe(df,name):
	grouped = df.groupby(by=[name])
	return grouped

def top_n(grouped,stat,n):
	for name,group in grouped:
		print str(group_name) + ": " + str(name) + "\n"
		print group.sort_values(stat).head(n)
		print "\n\n"

path = options.input
group_name = options.group_by
stat = options.stat
n = options.n
plot_params = [22,16,16,10,12]
a = 'overlaps'
b = 'combo_score'

df = parse_in(path)
df = rename_labels(df)
grouped = group_dataframe(df,group_name)
top_n(grouped,stat,n)
#make_scatter(df,group_name,a,b,plot_params)
#make_barplot(df,group_name,stat,plot_params)
