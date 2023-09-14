#!/usr/bin/env python2
#This script makes dihedral graphs for alpha,..,epsilon
 
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches  as patches
from   matplotlib import colors, ticker, cm, rc
from   matplotlib.colors import ListedColormap
from   scipy.interpolate import spline
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   mpl_toolkits.axes_grid1 import ImageGrid
import sys
import glob
import time
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from numpy.random import randn
from numpy import genfromtxt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as tck


csfont = {'fontname':'Comic Sans MS'}
hfont = {'fontname':'fantasy'}


font0 = FontProperties()
families = ['cursive']
font1 = font0.copy()
font1.set_size('large')

#inforna BBB
data_1 = open(sys.argv[1], 'r')
lines_1 = data_1.readlines()[1:-1]
data_1.close()

#drugbank new lib
#data_2 = open(sys.argv[2], 'r')
#lines_2 = data_2.readlines()[1:-1]
#data_2.close()

#drugbank ROBIN
#data_3 = open(sys.argv[3], 'r')
#lines_3 = data_3.readlines()[1:-1]
#data_3.close()

data1 = []
#data2 = []
#data3 = []

for line in lines_1:
    p = line.split()
    if float(p[0]) >=0:
        data1.append(float(p[0]))

#for line in lines_2:
#    p = line.split()
    #print('line', line)
#    data2.append(float(p[0]))

#for line in lines_3:
#    p = line.split()
    #print('line', line)
#    data3.append(float(p[0]))

xv1 = np.array(data1)
xv1_mean = np.mean(xv1)
xv1_std = np.std(xv1)


#yv1 = np.array(data2)
#yv1_mean = np.mean(yv1)
#yv1_std = np.std(yv1)

#zv1 = np.array(data3)
#zv1_mean = np.mean(zv1)
#zv1_std = np.std(zv1)


# Create a figure instance
fig, ax = plt.subplots()

# Create an axes instance
xticklabels = ['MW','MW','MW']
ax.set_xticks([1,2,3])
ax.set_xticklabels(xticklabels)

# Create the boxplot
#ax.violinplot(data_to_plot,showmedians=True)
sns.violinplot(data=[xv1],inner="quart", linewidth=1,cut=1,showextrema=False)

ax.set_xticklabels(['Inforna'],fontsize=15)


# Hide the right and top spines
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)

ax.set_ylabel('BBB',fontsize=15)
#ax.set_xticks([])
#ax.tick_params(axis='y', labelsize=14)
#ax.yaxis.labelpad = 10
#ax.set_xticks(ind)
#ax.set_xticklabels()

#fig.tight_layout()
plt.savefig('BBB.png',dpi=300)
plt.show()




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
