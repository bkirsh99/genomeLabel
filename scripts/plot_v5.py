#!/usr/bin/python
import glob
import sys
import math
from optparse import OptionParser
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
  
parser = OptionParser()

parser.add_option("--iq",
                  #"--input_query_dir",
                  dest="input_dir",
                  help="Input query directory")

parser.add_option("--ib",
				  #"--input_background_dir
				  dest="background_dir",
				  help="Background directory")
				  
parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-m",
				  "--map",
				  dest="map_file",
				  help="ENCODE url to factor file name")

parser.add_option("--stat",
                  dest="stat",
                  default="combo",
                  help="Stat to plot (odds, sig, combo) (Default: combo)")

(options, args) = parser.parse_args()

map = {}

if options.stat not in ['odds', 'sig', 'combo']:
    parser.error('Stat "' + options.stat + '" not supported')

f1 = open(options.map_file, 'r')
for line in f1:
	A = line.rstrip().split()
	match = A[1].split('/')[1].split('.')[0]
	map[match] = A[0]
	
def toMatrix(dir,c):
	results ={}
	X_names = []
	Y_names = []
	
	for file_name in glob.glob(dir):
		if os.stat(file_name).st_size != 0:
			find = file_name.split('/')[-1].split('.')[0]
			factor_name = ''
			for match in map.keys():
				if find == match:
					factor_name = map[match]
			Y_names.append(factor_name)
			results[factor_name] = {}
			f = open(file_name, 'r')
			next(f) 
			for l in f:
				A = l.rstrip().split('\t')
				odds = float(A[3]) + c
				sig = float(A[4]) + c
				combo = float(A[7]) + c
				gene_name = A[0].split('/')[-1].split('.')[:-1][0]
				X_names.append(gene_name) if gene_name not in X_names else X_names

				if options.stat == 'sig':
					results[factor_name][gene_name] = sig
				elif options.stat == 'combo':
					results[factor_name][gene_name] = combo
				elif options.stat == 'odds':
					results[factor_name][gene_name] = odds

	D = np.zeros([len(Y_names),len(X_names)])
	Y_names = sorted(Y_names, key=str.lower)
	X_names = sorted(X_names, key=str.lower)

	y = 0
	for yval in Y_names:
		x = 0
		for xval in X_names:
			D[y,x] = results[yval][xval]
			x+=1
		y+=1
	return (D,X_names,Y_names)

const = 0.001
(A,X_names,Y_names) = toMatrix(options.input_dir,const)
B = toMatrix(options.background_dir,const)[0]
#print A
#print B
ratio = np.divide(A,B)
#B_inv = np.array([1/i if i!=0 else 0 for i in B], dtype="float")
#data = np.transpose(A*B_inv)
data = np.transpose(ratio)
#print ratio
#print data
row_labels = X_names
column_labels = Y_names
row_size = 10
#column_size = 30
column_size = 10

fig =  matplotlib.pyplot.figure(figsize=(column_size,row_size), \
                                    dpi=300, \
                                    facecolor='black')

ax = fig.add_subplot(1,1,1,facecolor='k')

from matplotlib import colors as mcolors
from matplotlib.colors import Normalize

_seismic_data = ( (0.0, 0.0, 0.3),
                  (0.0, 0.0, 1.0),

                  (1.0, 1.0, 1.0),

                  (1.0, 0.0, 0.0),
                  (0.5, 0.0, 0.0))


hm = mcolors.LinearSegmentedColormap.from_list( \
        name='red_white_blue', \
        colors=_seismic_data, N=256)

class MidpointNormalize(Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		Normalize.__init__(self, vmin, vmax, clip)
	def __call__(self, value, clip=None):
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))

from pylab import get_cmap

#print data.min(), data.max()
plt.pcolor(data,
           #cmap=get_cmap("Reds"))
           #cmap=get_cmap("Blues"))
           cmap=hm,
           norm = MidpointNormalize(midpoint=0))

ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
plt.xticks(np.arange(0.5,len(column_labels)+0.5,1.0),column_labels,rotation=90, fontsize=10, color='white')
plt.yticks(np.arange(0.5,len(row_labels)+0.5,1.0),row_labels,fontsize=10, color='white')
plt.ylim((0,len(row_labels)))
plt.xlim((0,len(column_labels)))
cbar = plt.colorbar(fraction=0.046, pad=0.04, ticks=[data.min(), 0, data.max()])
cbar.ax.tick_params(labelsize=10, color='white', labelcolor='white')
plt.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
            transparent=True)
