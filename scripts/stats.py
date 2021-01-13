#!/usr/bin/python
#!pip install brewer2mpl

import glob
import sys
import math
from optparse import OptionParser
import os
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from bioinfokit import analys, visuz
import openpyxl

parser = OptionParser()

parser.add_option("-i",
                  "--input",
                  dest="input",
                  help="Input results directory or file")

parser.add_option("-o",
                  "--output",
                  dest="output",
                  help="Output file name")

parser.add_option("--do",
				  dest="do",
				  default="summary",
				  help="Output to create (summary, clustermap)")

parser.add_option("--stat",
                  dest="stat",
                  default="combo_score",
                  help="Stat to describe (query_file_size, overlaps, odds_ratio, fishers_two_tail, fishers_left_tail, fishers_right_tail, combo_score) (Default: combo_score)")

parser.add_option("--highlight",
 				   dest="highlight",
 				   type="str",
 				   help="Labels to highlight")
   
parser.add_option("--x_stat",
				  dest="x_stat",
				  help="X-axis stat")

parser.add_option("--y_stat",
                 dest="y_stat",
				 help="Y-axis stat")

parser.add_option("--x_size",
				  dest="x_size",
                  type="int",
                  #default=30,
                  help="Figure x size (Default 10)")

parser.add_option("--y_size",
                  dest="y_size",
                  type="int",
                  #default=20,
                  help="Figure x size (Default 30)")

parser.add_option("--by",
				  dest="group_by",
                  default="query",
                  help="Output by query or index annotations")

parser.add_option("--labels",
                  dest="labels",
                  default=True,
                  help="Label x and y axis");

parser.add_option("--xlabels",
                  dest="xlabels",
                  default=True,
                  help="Label x axis");

parser.add_option("--ylabels",
                  dest="ylabels",
                  default=True,
                  help="Label y axis");

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

valid_stats = ['query_file_size', 'overlaps', 'odds_ratio', 'fishers_two_tail', 'fishers_left_tail', 'fishers_right_tail', 'combo_score']
if options.stat not in valid_stats:
    parser.error('Stat "' + options.stat + '" not supported')

def normalize(df):
	result = df.copy()
	for feature_name in df.columns:
		max_value = df[feature_name].max()
		min_value = df[feature_name].min()
		result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
	return result

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

def make_heatmap(arr,stat,highlight,dim,out,x,y):
	#data = to_2darray(df,stat)
	#print(data)
	#data_norm = normalize(data)
	#print(data_norm)
	#data = data_norm.dropna()
	#visuz.gene_exp.hmap(df=data, dim=dim, tickfont=(6, 4),figname=out,xlabel=x,ylabel=y,zscore=0)
	print(arr)
	if stat in ['fishers_two_tail', 'fishers_left_tail', 'fishers_right_tail']:
		arr = np.negative(np.log10(arr))
	print(arr)
	#mask_na = 0.000666
	#arr = arr.fillna(mask_na)
	#arr = arr.add(0.01)
	#arr = arr.loc[(arr!=0).any(axis=0)]
	arr = arr.replace(0, np.nan)
	arr = arr.dropna(how='all', axis=0)
	arr = arr.replace(np.nan, 0)
	print(arr)
	color_dict = {}
	palette = sns.color_palette()
	for col in arr.columns:
		if col in highlight:
			color_dict[col] = palette[0]
		else:
			color_dict[col] = palette[1]
	color_rows = pd.Series(color_dict)
	#row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (arr.values, arr.values.T))
	#print(col)
	#arr = normalize(arr)
	#arr = arr.dropna()
	#arr[~arr.isin([np.nan, np.inf, -np.inf]).any(1)]
	#arr = arr.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
	#sns.clustermap(arr,figsize=dim,cmap='coolwarm',col_colors=[color_rows],standard_scale=0,row_linkage=row_linkage,col_linkage=col_linkage)
	sns.clustermap(arr,figsize=dim,cmap='coolwarm',col_colors=[color_rows],standard_scale=0)
	plt.savefig(out)
	
def parse_in(path):
	results = ''
	ith_file = 0
	files = []

	if os.path.isdir(path):
		for file in glob.glob(os.path.join(path, '*')):
			(results, ith_file, files) = filter_cluster(file,results,ith_file,files)
	elif os.path.isfile(path):
		(results, ith_file, files) = filter_cluster(path,results,ith_file,files)

	res_table = pd.read_csv(io.StringIO(results),sep="\t",header=None)
	res_table.columns = ['query', 'index', 'query_file_size', 'overlaps', 'odds_ratio', 'fishers_two_tail', 'fishers_left_tail', 'fishers_right_tail', 'combo_score']
	return rename_labels(res_table)

def filter_cluster(file,str,n,arr):
	if os.stat(file).st_size != 0:
		with open(file, 'r') as f:
			result = f.read().replace('\t\n', '\n')
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
	df['index'] = df['index'].apply(lambda x: os.path.basename(x))
	df = df.replace(regex=r'.bed.gz.*$',value='')
	return df

def to_2darray(table,val):
	#df = table.set_index('index')
	df = table.pivot_table(index='index', columns='query', values=val) 
	return df

def group_dataframe(df,name):
	grouped = df.groupby(by=[name])
	return grouped

def top_n(grouped,stat,n):
	for name,group in grouped:
		print(str(group_name) + ": " + str(name) + "\n")
		print(group.sort_values(stat).head(n))
		print("\n\n")

def print_summary(table,val,out):
	writer = pd.ExcelWriter(out,engine='xlsxwriter')
	df = table.set_index('index')
	for value in valid_stats:
		pivot = df.pivot(columns='query',values=value)
		pivot.to_excel(writer,sheet_name=value)
		pivot.describe().to_excel(writer,sheet_name='summary_' + value)
	writer.save()
	return df

path_in = options.input
path_out = options.output 
group_name = options.group_by
stat = options.stat
n = options.n
#width = options.x_size if options.x_size else 5+0 if len(data.columns)<50 else (len(data.columns)-50)/100
#row_cutoff = 1000
#height = options.y_size if options.y_size else 15+0 if len(data)<row_cutoff else (len(data)-row_cutoff)/75.0
#dim = (width,height)
dim = (options.x_size,options.y_size)
xlabels = options.xlabels
ylabels = options.ylabels
highlight = options.highlight.split(",")
print(highlight)
plot_params = [22,16,16,10,12]
a = 'overlaps'
b = 'combo_score'

tbl = parse_in(path_in)
#df = rename_labels(df)
#print(tbl)
arr = to_2darray(tbl,stat)
#print(arr)
#summ = print_summary(tbl,stat,path_out)
#print(summ)
#grouped = group_dataframe(df,group_name)
#top_n(grouped,stat,n)
#width = options.x_size if options.x_size else 5+0 if len(arr.columns)<50 else (len(arr.columns)-50)/100
#row_cutoff = 1000
#height = options.y_size if options.y_size else 15+0 if len(arr)<row_cutoff else (len(arr)-row_cutoff)/75.0
#dim = (width,height)
#print(dim)
#print(len(arr.columns))
#print(len(arr))

#make_scatter(df,group_name,a,b,plot_params)
#make_barplot(df,group_name,stat,plot_params)
make_heatmap(arr,stat,highlight,dim,path_out,xlabels,ylabels)
