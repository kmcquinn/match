# -*- coding: utf-8 -*-
"""
Created on Sun May 01 20:28:46 2016

@author: Anna

Directory structure: assumes galaxy directory is the working directory, input_data/ should contain pars, phot and fake files
This code generates a CMD plot in the same directory where the program, CMDscript.py, is
Use optional -phot, -fake flags to specify their location, or change the default setting in parse_options()
Assumes the default values of the column numbers for fake data to be 2,4,9 and 11; specify the values in the bash script or edit the default values
Give any name for .png file after -fake=... The file gets saved as [filename].png
Syntax: python CMDscript.py [path to working directory, ending in /] -phot=[path from working directory to phot file] -fake=...

e.g. python CMDscript.py /work/04316/kmcquinn/wrangler/shield/galaxies/agc223254/ -phot=phot/a223254_tilted_ellipse.phot3 2 4 9 11 -fake=input_data/fake AGC223254
The output files (.png files) are stored in the same directory as this program's

"""
import matplotlib
matplotlib.use("Agg")# Comment this line if not running on tacc
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import pylab as pl
import sys

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

import argparse


#def parse_options():
    # Creates argument parser
parser = argparse.ArgumentParser(description='Process input parameters.')
    # Defines required arguments used on the command line
parser.add_argument('galaxydir',action='store',help='path to galaxy directory (note: this is the whole path ending in the galaxy directory)')
parser.add_argument('-phot',dest='phot',action='store',help="location of phot file from galaxy directory",default='input_data/phot')
parser.add_argument('col1',action='store',help="column number for V",default='2', type=int)
parser.add_argument('col2',action='store',help="column number for Verr",default='4', type = int)
parser.add_argument('col3',action='store',help="column number for I",default='9',type = int)
parser.add_argument('col4',action='store',help="column number for Ierr",default='11',type=int)
parser.add_argument('-fake',dest='fake',action='store',help="location of fake file from galaxy directory",default='input_data/fake')
parser.add_argument('gname',action = 'store', help="names the galaxy", type=str)

    # Parses through the arguments and saves them within the keyword args
args = parser.parse_args()
    
#plot_title = sys.argv[1]
#redfilter = sys.argv[2]
#bluefilter = sys.argv[3]
#plot_title = "AGC238890"
 # Parses through command line arguments
#args = parse_options()
galdir = args.galaxydir
#col_1 = args.col1
#col_2 = args.col2
#col_3 = args.col3
#col_4 = args.col4

    # Defines the location of phot, and fake

phot = galdir + args.phot
fake = galdir + args.fake
    #print filtername(pars)
gal_name = args.gname    
    
plot_title= gal_name
redfilter = "F814W"
bluefilter = "F606W"
print(str(galdir))

#Real Data
"""Parsers the path to phot file"""

dat = np.genfromtxt(phot)
#dat = np.genfromtxt(args.phot+'')

#remove all NaNs
dat = dat[~np.isnan(dat).any(axis=1)]

"""Parsers the column numbers for fake"""
#V = np.array(dat[:,2])
V= np.array(dat[:,args.col1])
#Verr = np.array(dat[:,4])
Verr = np.array(dat[:,args.col2])
#I = np.array(dat[:,9])
I = np.array(dat[:,args.col3])
#Ierr = np.array(dat[:,11])
Ierr = np.array(dat[:,args.col4])
VmI = V-I

VmIerr = (((Verr)**2)+(((Ierr)**2)))**.5
print VmIerr

#fake error stuff

fdat = np.loadtxt(fake)
#fdat = np.loadtxt(args.fake)

fdat = np.asarray([d for d in fdat if not 99.999 in d])
fVerr = np.array(fdat[:,2])
fIerr = np.array(fdat[:,3])
fI = np.array(fdat[:,1])
fV = np.array(dat[:,0])
fVmIerr = (fVerr**2 + fIerr**2)**0.5
print fVmIerr

#Here I am finding the max and min values of the data
#We will use this to automate the figsize thing 
maxV = np.amax(V) + .3
print max(V)
maxI = np.amax(I) + .3
minV = np.amin(V) + 1.3
minI = np.amin(I) + 1.3
meanVmI = np.mean(VmI)
maxVmI = (meanVmI) + 5*np.std(VmI)
minVmI = (meanVmI) - 4*np.std(VmI)
Ierrup = np.around(maxI - .5)
Verrup = np.around(maxV - .5)
Ierrlow = np.around(minI + .5)
Verrlow = np.around(minV + .5)
errx = maxVmI -.3
print maxVmI
print errx
#maxVmI = np.amax(VmI) + .2
#minVmI = np.amin(VmI) - .2
#print minI
#print min(I)



#Making plain scatter plot
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1) 
ax.scatter(VmI, I, color = 'black', s = 2)
ax.set_ylim(maxI, minI)
ax.set_xlim(minVmI, maxVmI)
ax.set_xlabel(bluefilter + '-' + redfilter + ' (mag)')
ax.set_ylabel(redfilter + ' (mag)')
fig.text(.18, .85, 'Number of stars: '+str(len(V)))

majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = AutoMinorLocator(10)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)

plt.suptitle(plot_title)
#ax.set_title(plot_title)

#adding error to the plot
#c = 3
c = errx
#location of the line
errlist = []
ylist = []
Verrlist = []
VmIerrlist = []
#Ilist = range(18, 27)
Ilist = range(int(Ierrlow), int(Ierrup))

for a in Ilist:
    Iw = np.where(I > a)
    Iwr = np.where(I[Iw] < a+1)
    fIw = np.where(fI > a)
    fIwr = np.where(fI[fIw] < a+1)
    Ierravg = ((np.mean(Ierr[Iwr]))**2 + (np.mean(abs(fIerr[fIwr]))**2))**.5
    errlist.append(Ierravg)
    ylist.append(a+.5)
    VmIerravg = (np.mean(VmIerr[Iwr])**2 + np.mean(abs(fVmIerr[fIwr]))**2)**.5
    VmIerrlist.append(VmIerravg)
    
xlist = c*np.ones_like(ylist)
plt.errorbar(xlist, ylist, xerr = VmIerrlist, yerr= errlist, fmt = '.', capsize=0)


#adding in the contour script
def multidigitize(VmI,I,binsVmI,binsV):
    dVmI = np.digitize(VmI.flat, binsVmI)
    dI = np.digitize(I.flat, binsV)
    return dVmI,dI

def linlogspace(xmin,xmax,n):
    return np.logspace(np.log10(xmin),np.log10(xmax),n)

#here's the contour actual values
def adaptive_param_plot(VmI,I,
                        bins=3,
                        threshold=2,
                        marker='.',
                        marker_color=None,
                        ncontours=5,
                        fill=False,
                        mesh=False,
                        contourspacing=linlogspace,
                        mesh_alpha=0.5,
                        norm=None,
                        axis=None,
                        cmap=None,
                        **kwargs):
    if axis is None:
        axis = pl.gca()
        axis.set_ylim(28, 18)
    ok = np.isfinite(VmI)*np.isfinite(I)
    
    if hasattr(bins, 'ndim') and bins.ndim == 2:
        nbinsVmI, nbinsI = bins.shape[0]-1, bins.shape[1]-1
    else:
        try:
            nbinsVmI = nbinsI = len(bins)-1
        except TypeError:
            nbinsVmI = nbinsI = bins
    H, bVmI, bI = np.histogram2d(VmI[ok], I[ok], bins = bins)
    
    dVmI, dI = multidigitize(VmI[ok], I[ok], bVmI, bI)
    
    plottable = np.ones([nbinsVmI+2, nbinsI+2], dtype = 'bool')
    plottable_hist = plottable[1:-1, 1:-1]
    assert H.shape == plottable_hist.shape
    
    plottable_hist[H > threshold] = False
    
    H[plottable_hist] = 0
    
    toplot = plottable[dVmI, dI]
    
    cVmI = (bVmI[1:]+bVmI[:-1])/2
    cI = (bI[1:]+bI[:-1])/2
    levels = contourspacing(threshold-0.5, H.max(), ncontours)
    
    if cmap is None:
        cmap = plt.cm.get_cmap()
        cmap.set_under((0,0,0,0))
        cmap.set_bad((0,0,0,0))
    
    if fill:
        con = axis.contourf(cVmI, cI, H.T, levels= levels, norm = norm, cmap = cmap,  **kwargs)
    else: 
        con = axis.contour(cVmI, cI, H.T,levels=levels,norm=norm,cmap=cmap,**kwargs) 
    if mesh: 
        mesh = axis.pcolormesh(bVmI, bI, H.T, **kwargs)
        mesh.set_alpha(mesh_alpha)
        #Is there a way to add lines w the contour levels?
    
    if 'linestyle' in kwargs:
        kwargs.pop('linestyle')
        
    #if i wanted to plot the scatter from this script intstead, but I can't make it look as nice
   # axis.plot(VmI[ok][toplot],
    #          I[ok][toplot],
     #         linestyle='none',
      #        marker=marker,
       #       markerfacecolor=marker_color,
        #      markeredgecolor=marker_color,
         #     **kwargs)
    
    return cVmI, cI, H, VmI[ok][toplot], I[ok][toplot]

adaptive_param_plot(VmI, I, bins = 100, fill = True, ncontours = 7, threshold = 2, axis = ax)





#SECOND PLOT




#Making plain scatter plot
#fig = plt.figure()

ax1 = fig.add_subplot(1, 2, 2) 
ax1.scatter(VmI, V, color = 'black', s = 2)
ax1.set_ylim(maxV, minV)
ax1.set_xlim(minVmI, maxVmI)
#ax.set_ylim(29, 18)
#ax.set_xlim(-.98, 3.48)
ax1.set_xlabel(redfilter + '-' + bluefilter + ' (mag)')
ax1.set_ylabel(bluefilter + ' (mag)')
ax1.yaxis.set_label_coords(-0.12, 0.5)
#fig.text(.6, .85, 'Number of stars: '+str(len(V)))
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = AutoMinorLocator(10)
ax1.yaxis.set_major_locator(majorLocator)
ax1.yaxis.set_major_formatter(majorFormatter)
ax1.yaxis.set_minor_locator(minorLocator)


ax.get_shared_y_axes().join(ax, ax1)
ax.set_yticklabels([])
ax1.autoscale()

#ax.set_title(plot_title)


#adding error to the plot
#c = 3
c = errx
#location of the line
errlist = []
ylist = []
Verrlist = []
VmIerrlist = []
#Ilist = range(18, 28)
Ilist = range(int(Verrlow), int(Verrup))

for a in Ilist:
    Iw = np.where(I > a)
    Iwr = np.where(I[Iw] < a+1)
    fIw = np.where(fI > a)
    fIwr = np.where(fI[fIw] < a+1)
    Ierravg = ((np.mean(Ierr[Iwr]))**2 + (np.mean(abs(fIerr[fIwr])))**2)**.5
    errlist.append(Ierravg)
    ylist.append(a+.5)
    VmIerravg = (np.mean(VmIerr[Iwr])**2 + np.mean(fVmIerr[fIwr])**2)**.5
    VmIerrlist.append(VmIerravg)

xlist = c*np.ones_like(ylist)
plt.errorbar(xlist, ylist, xerr = VmIerrlist, yerr= errlist, fmt = '.', capsize = 0)


#adding in the contour script
def multidigitize(VmI,V,binsVmI,binsV):
    dVmI = np.digitize(VmI.flat, binsVmI)
    dV = np.digitize(V.flat, binsV)
    return dVmI,dV

def linlogspace(xmin,xmax,n):
    return np.logspace(np.log10(xmin),np.log10(xmax),n)

#here's the contour actual values
def adaptive_param_plot(VmI,V,
                        bins=3,
                        threshold=2,
                        marker='.',
                        marker_color=None,
                        ncontours=5,
                        fill=False,
                        mesh=False,
                        contourspacing=linlogspace,
                        mesh_alpha=0.5,
                        norm=None,
                        axis=None,
                        cmap=None,
                        **kwargs):
    if axis is None:
        axis = pl.gca()
        axis.set_ylim(28, 18)
    ok = np.isfinite(VmI)*np.isfinite(V)
    
    if hasattr(bins, 'ndim') and bins.ndim == 2:
        nbinsVmI, nbinsV = bins.shape[0]-1, bins.shape[1]-1
    else:
        try:
            nbinsVmI = nbinsV = len(bins)-1
        except TypeError:
            nbinsVmI = nbinsV = bins
    H, bVmI, bV = np.histogram2d(VmI[ok], V[ok], bins = bins)
    
    dVmI, dV = multidigitize(VmI[ok], V[ok], bVmI, bV)
    
    plottable = np.ones([nbinsVmI+2, nbinsV+2], dtype = 'bool')
    plottable_hist = plottable[1:-1, 1:-1]
    assert H.shape == plottable_hist.shape
    
    plottable_hist[H > threshold] = False
    
    H[plottable_hist] = 0
    
    toplot = plottable[dVmI, dV]
    
    cVmI = (bVmI[1:]+bVmI[:-1])/2
    cV = (bV[1:]+bV[:-1])/2
    levels = contourspacing(threshold-0.5, H.max(), ncontours)
    
    if cmap is None:
        cmap = plt.cm.get_cmap()
        cmap.set_under((0,0,0,0))
        cmap.set_bad((0,0,0,0))
    
    if fill:
        con = axis.contourf(cVmI, cV, H.T, levels= levels, norm = norm, cmap = cmap,  **kwargs)
    else: 
        con = axis.contour(cVmI, cV, H.T,levels=levels,norm=norm,cmap=cmap,**kwargs) 
    if mesh: 
        mesh = axis.pcolormesh(bVmI, bV, H.T, **kwargs)
        mesh.set_alpha(mesh_alpha)
        #Is there a way to add lines w the contour levels?
    
    if 'linestyle' in kwargs:
        kwargs.pop('linestyle')

    
    return cVmI, cV, H, VmI[ok][toplot], I[ok][toplot]

adaptive_param_plot(VmI, V, bins = 100, fill = True, ncontours = 7, threshold = 2, axis = ax)

#plt.savefig('CMD.png')
plt.savefig(gal_name+'.png')




