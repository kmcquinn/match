#original code Hess_plot.py, revised by Anthony Pahl 11/30/2016
#
#syntax: CMD_res_SF_Z_plot.py out2.cmd out2.final output2.png
#
#plots observed CMD, modelled CMD, and residual based on calcsfh
#outputs. also plots SF, SFR, and Z vs. time for said fits.
#attempts to recreate 'plot_cmd' output images from match.

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import sys

cmdfile_name = sys.argv[1]
outfile_name = sys.argv[2]
plot_name = sys.argv[3]
figsize=(8, 8)
subsize = 0.25


subx1, suby1 = 0.05, 0.05
subx2, suby2 = 0.38, 0.38


subx3, suby3 = 0.72, 0.72
interpolation = 'hanning'
xlabel, ylabel = 'F606W - F814W', 'F606W'

def makecmap(arr):
    x = np.linspace(arr.min(), arr.max(), 50)
    steps = (x - x.min()) / (x.max() - x.min())
    numerator = np.arcsinh(x) - np.arcsinh(x.min())
    denominator = np.arcsinh(x.max()) - np.arcsinh(x.min())
    mapping = 1 - numerator / denominator
    cdict = {}
    for key in ('red', 'green', 'blue'):
        cdict[key] = np.vstack([steps, mapping, mapping]).transpose()
    return cl.LinearSegmentedColormap('new_colormap', cdict, N=1024)

def plot_HessD(fig, arr, subx1, suby2, subsize, extent, interpolation,
               title, xlabel, ylabel,col_bool):
    if col_bool:
        cm = matplotlib.cm.jet#makecmap(arr)
        nm = None
        cmap_ofs=0.
    else:
        cm=makecmap(arr)
        nm=None
        cmap_ofs=0.05
    #cm.set_gamma(0.2)
    ax = fig.add_axes([subx1, suby2, subsize, subsize])
    ax.imshow(arr, cmap=cm, aspect='auto', extent=extent,
              interpolation=interpolation, norm=nm)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xticklabels(ax.get_xticks(), size=8)
    ax.set_yticklabels(ax.get_yticks(), size=8)
    #cax = fig.add_axes([subx1+0.13, suby2+0.35, subsize-0.15, subsize-0.365])
    cax = fig.add_axes([subx1+0.07+cmap_ofs, suby2+0.225,
                        subsize-0.10-cmap_ofs, subsize-0.235])
    cb_arr = np.linspace(arr.min(), arr.max(), 1e2).reshape(1, 1e2)
    cax_limits = [np.floor(arr.min()), np.ceil(arr.max()), 0, 1]
    cax.imshow(cb_arr, cmap=cm, aspect='auto', extent=cax_limits, interpolation='nearest', norm=nm)
    cax.plot([arr.min(), arr.min()], [0, 1], 'k-')
    cax.plot([arr.max(), arr.max()], [0, 1], 'w-')
    cax.axis(cax_limits)
    cax.yaxis.set_visible(False)
    cax.xaxis.tick_bottom()
    cax.xaxis.set_major_locator(plt.LinearLocator(5))
    cax.set_xticklabels(cax.get_xticks(), size=5)
    #cax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    return ax, cax

def plot_Scatter(fig, x, y, subx1, suby2, subsize, title, xlabel, ylabel,linestyle, errlow=[], errhigh=[], bound=False):
	ax = fig.add_axes([subx1, suby2, subsize, subsize])
	if len(errlow) == 0:
		ax.plot(x, y,'k', linestyle = linestyle)
	else:
		ax.errorbar(x, y, yerr=[errlow, errhigh], fmt='k', capsize=0, ecolor='k', elinewidth=0.7)
	ax.set_title(title, fontsize=8)
	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	if bound == True:
		ax.set_xlim(14.0, min(x))
	else:
		ax.set_xlim(max(x),min(x))
	ax.set_xticklabels(ax.get_xticks(), size=8)
	ax.set_yticklabels(ax.get_yticks(), size=8)
	return ax
    
def plot_Hist(fig, x, y, subx1, suby2, subsize, title, xlabel, ylabel,linestyle, errlow=[], errhigh=[], bound=False):
	#x: array with all x values
	#y: array with corrosponding y values
	ax = fig.add_axes([subx1, suby2, subsize, subsize])
	ax.step(x, y, color='k')
	ax.errorbar(x,y,yerr=[errlow,errhigh], fmt='none', capsize=0, ecolor='k', elinewidth=0.7)
	ax.set_title(title, fontsize=8)
	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	if bound == True:
		ax.set_xlim(14.0, min(x))
	else:
		ax.set_xlim(max(x),min(x))
	ax.set_xticklabels(ax.get_xticks(), size=8)
	ax.set_yticklabels(ax.get_yticks(), size=8)
	return ax

def plot_fill(fig, x, y, subx1, suby2, subsize, title, xlabel, ylabel,linestyle, errlow=[], errhigh=[]):
    ax = fig.add_axes([subx1, suby2, subsize, subsize])
    ax.plot(x, y,'k', linestyle = linestyle)
    ax.fill_between(x, errlow, errhigh, facecolor='k', alpha=0.5)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlim(max(x),min(x))
    ax.set_xticklabels(ax.get_xticks(), size=8)
    ax.set_yticklabels(ax.get_yticks(), size=8)
    return ax

cmdfile = open(cmdfile_name, 'r')
line = cmdfile.readline()
shape = [int(i) for i in cmdfile.readline().split()]
shape = shape[1:3]
line_col, line_mag = cmdfile.readline(), cmdfile.readline()
magbins, colbins = np.zeros(shape[0]), np.zeros(shape[1])
obs_arr, mod_arr = np.zeros(shape), np.zeros(shape)
res_arr, sig_arr = np.zeros(shape), np.zeros(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        line = cmdfile.readline().split()
        mag, col, obs, mod, res, sig = [float(n) for n in line[:6]]
        colbins[j] = col
        obs_arr[i, j], mod_arr[i, j] = obs, mod
        res_arr[i, j], sig_arr[i, j] = res, sig
    magbins[i] = mag
cmdfile.close()

outfile = open(outfile_name, 'r')
lage_arr=np.zeros(71)
sfr_arr=np.zeros(71)
sfr_low=np.zeros(71)
sfr_high=np.zeros(71)
csfr_arr=np.zeros(71)
To_arr = np.zeros(71)
Tf_arr = np.zeros(71)
met_arr = np.zeros(71)
met_up = np.zeros(71)
met_low = np.zeros(71)
csf_arr = np.zeros(71)
csf_up = np.zeros(71)
csf_low = np.zeros(71)
line =outfile.readline()
totalSF = float(line.split()[1])
for i in range(71):
    line = outfile.readline().split()
    age, Tf, x, sfr, y, z, met, metup, metdown, a, b, c, csf, csfup, csfdown = [float(n) for n in line[:15]]
    lage_arr[i]=age
    sfr_arr[i]=sfr
    sfr_high[i] = y
    sfr_low[i] = z
    met_arr[i]=met
    met_up[i] = metup
    met_low[i] = metdown
    csf_arr[i] = csf
    csf_up[i] = csfup
    csf_low[i] = csfdown
outfile.close

age_arr = [10**(n)/1e9 for n in lage_arr]
for i in range(71):
    if i==0:
        csfr_arr[i]=0.
    else:
        csfr_arr[i]=csfr_arr[i-1]-sfr_arr[i]*(10**(lage_arr[i])-10**(lage_arr[i]-0.05))
#csfr_arr = csfr_arr + totalSF
csfr_arr = [(n+totalSF) / totalSF for n in csfr_arr]



xlabel = ' - '.join(line_col.replace('WFC','F').split('-'))
ylabel = line_mag.replace('WFC','F')

fig = plt.figure(figsize=figsize)

extent = [colbins[0], colbins[-1], magbins[-1], magbins[0]]
ax_obs, cax_obs = plot_HessD(fig, obs_arr, subx1, suby3, subsize, extent, interpolation, '(a) Observed CMD', xlabel, ylabel,1)

ax_mod, cax_mod = plot_HessD(fig, mod_arr, subx1, suby2, subsize, extent, interpolation, '(d) Modeled CMD', xlabel, ylabel,1)

#ax_res, cax_res = plot_HessD(fig, res_arr, subx1, suby1, subsize, extent, interpolation, '(c) Residual (Obs. - Mod.)', xlabel, ylabel)

ax_sig, cax_sig = plot_HessD(fig, sig_arr, subx1, suby1, subsize, extent, interpolation, '(g) Residual Significance', xlabel, ylabel,0)

ax_2 = plot_Hist(fig, lage_arr, sfr_arr, subx2, suby2, subsize, '(e) Star Formation Rate', 'log(age)','SFR','steps', sfr_low, sfr_high)

ax5 = plot_Hist(fig, age_arr, sfr_arr, subx3, suby2, subsize, '(f) Star Formation Rate', 'age (Gyr)','SFR','steps', sfr_low, sfr_high, True)

ax_1 = plot_Scatter(fig, lage_arr, csf_arr, subx2, suby3, subsize, '(b) Cumulative Star Formation', 'log(age)','Fraction of Stellar Mass','solid', csf_low, csf_up)

ax_4 = plot_Scatter(fig, age_arr, csf_arr, subx3, suby3, subsize, '(c) Cumulative Star Formation', 'age (Gyr)', 'Fraction of Stellar Mass','solid', csf_low, csf_up, True)


ax_3 = plot_Scatter(fig, lage_arr, met_arr, subx2, suby1, subsize, '(h) Metallicity', 'log(age)','Z','solid', met_low, met_up)

ax_6 = plot_Scatter(fig, age_arr, met_arr, subx3, suby1, subsize, '(i) Metallicity', 'age (Gyr)', 'Z','solid', met_low, met_up, True)

fig.savefig(plot_name, dpi=300, bbox_inches='tight')
