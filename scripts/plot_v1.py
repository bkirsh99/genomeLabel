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

parser.add_option("-i",
                  "--input_dir",
                  dest="input_dir",
                  help="Input directory")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("--stat",
                  dest="stat",
                  default="combo",
                  help="Stat to plot (odds, sig, combo) (Default: combo)")

(options, args) = parser.parse_args()
results = {}
X_names = []
Y_names = []

if options.stat not in ['odds', 'sig', 'combo']:
    parser.error('Stat "' + options.stat + '" not supported')

for file_name in glob.glob(options.input_dir):
    gene_name = file_name.split('/')[-1].split('.')[:-1][0]
    Y_names.append(gene_name) if gene_name not in Y_names else Y_names
    results[gene_name] = {}
    f = open(file_name, 'r')
    next(f)
    for l in f:
    	A = l.rstrip().split('\t')
    	odds = float(A[3])
    	sig = float(A[4])
    	combo = float(A[7])
    	factor_name = A[0].split('/')[-1].split('_')[:-1][0]
    	X_names.append(factor_name) if factor_name not in X_names else X_names

    	if options.stat == 'sig':
    		results[gene_name][factor_name] = sig
    	elif options.stat == 'combo':
    		results[gene_name][factor_name] = combo
    	elif options.stat == 'odds':
    		results[gene_name][factor_name] = odds

D=np.zeros([len(Y_names),len(X_names)])

y = 0
for yval in Y_names:
	x = 0
	for xval in X_names:
		D[y,x] = results[yval][xval]
		x+=1
	y+=1

column_labels = Y_names
row_labels = X_names
data = D

def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):

	"""Create a heatmap from a numpy array and two lists of labels.
Parameters
----------
data
A 2D numpy array of shape (N, M).
row_labels
A list or array of length N with the labels for the rows.
col_labels
A list or array of length M with the labels for the columns.
ax
A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
not provided, use current axes or create a new one.  Optional.
cbar_kw
A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
cbarlabel
The label for the colorbar.  Optional.
**kwargs
All other arguments are forwarded to `imshow`."""
	if not ax:
		ax = plt.gca()
	im = ax.imshow(data, **kwargs)
	cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
	cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
	ax.set_xticks(np.arange(data.shape[1]))
	ax.set_yticks(np.arange(data.shape[0]))
	ax.set_xticklabels(col_labels)
	ax.set_yticklabels(row_labels)
	ax.tick_params(top=True, bottom=False,labeltop=True, labelbottom=False)
	plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",rotation_mode="anchor")
	for edge, spine in ax.spines.items():
		spine.set_visible(False)
	ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
	ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
	ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
	ax.tick_params(which="minor", bottom=False, left=False)
	return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",textcolors=["black", "white"],threshold=None, **textkw):

	"""A function to annotate a heatmap.
Parameters
----------
im
The AxesImage to be labeled.
data
Data used to annotate.  If None, the image's data is used.  Optional.
valfmt
The format of the annotations inside the heatmap.  This should either
use the string format method, e.g. "$ {x:.2f}", or be a
`matplotlib.ticker.Formatter`.  Optional.
textcolors
A list or array of two color specifications.  The first is used for
values below a threshold, the second for those above.  Optional.
threshold
Value in data units according to which the colors from textcolors are
applied.  If None (the default) uses the middle of the colormap as
separation.  Optional.
**kwargs
All other arguments are forwarded to each call to `text` used to create
the text labels."""

	if not isinstance(data, (list, np.ndarray)):
		data = im.get_array()
	if threshold is not None:
		threshold = im.norm(threshold)
	else:
		threshold = im.norm(data.max())/2.
	kw = dict(horizontalalignment="center", verticalalignment="center")
	kw.update(textkw)
	if isinstance(valfmt, str):
		valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
	texts = []
	for i in range(data.shape[0]):
		for j in range(data.shape[1]):
			kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
			text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
			texts.append(text)

	return texts

#fig, ax = plt.subplots()
#im, cbar = heatmap(data, Y_names, X_names, ax=ax,
#                   cmap="YlGn", cbarlabel="gene vs. factor")
#texts = annotate_heatmap(im, valfmt="{x:.1f} t")

fig =  matplotlib.pyplot.figure(figsize=(30,10), \
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
           #cmap=get_cmap("Rlues"))
           cmap=hm,
           norm = MidpointNormalize(midpoint=0))

ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
plt.xticks(np.arange(0.5,len(X_names)+0.5,1.0),X_names,rotation=90, fontsize=10, color='white')
plt.yticks(np.arange(0.5,len(Y_names)+0.5,1.0),Y_names,fontsize=10, color='white')
plt.ylim((0,len(Y_names)))
plt.xlim((0,len(X_names)))
cbar = plt.colorbar(fraction=0.046, pad=0.04, ticks=[data.min(), 0, data.max()])
cbar.ax.tick_params(labelsize=10, color='white', labelcolor='white')
plt.savefig(options.output_file,bbox_inches='tight',\
            facecolor=fig.get_facecolor(),\
            transparent=True)
