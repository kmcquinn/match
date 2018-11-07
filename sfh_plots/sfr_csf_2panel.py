"""
sfr_csf_2panel.py

Created by Zili Shen 2018

This script creates 2-panel histogram: 
left, SFR(t) with log time bins
right, CSF(t) with linear time axis and redshift axis on top

Input: out.final
Output: one png file called "name.png"

To use: modify galname to change plot title
modify name to change output file name
modify data to point to the desired out.final file
"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from astropy.cosmology import FlatLambdaCDM

galname = ('Star Formation History of KDG215')
data = np.genfromtxt('/Users/zilishen/core-cusp/12878_KDG215/out_c4p_1.final', skip_header = 1, unpack = True)
name = 'KDG215_sfh_two_panel'

fig = plt.figure(figsize=(16, 8))
outer = gridspec.GridSpec(1, 1, wspace=0.15, hspace=0.2)

ax_list=[]
for i in range(1):
    inner = gridspec.GridSpecFromSubplotSpec(1, 2,
                    subplot_spec=outer[i], wspace=0.2, hspace=0.1)
    ax_list.append([])

    for j in range(2):
        ax = plt.Subplot(fig, inner[j])
        ax_list[i].append(ax)
        fig.add_subplot(ax)
        
### Define plot
sfh_plot = ax_list[0][0]
lin_plot = ax_list[0][1]

### Data for linear plot
log_age = np.append(data[0],10.15)
lin_age = 10.**log_age / 1.0e9
lin_age[-1] = 13.7

csf = np.append(data[12],0)
csf_up = np.append(data[13],0)
csf_low = np.append(data[14],0)

### Log space time bins
t_begin = np.array(data[0])
t_end = np.array(data[1])
t_width = t_end - t_begin
t_mid = t_begin + 0.5*t_width

### Get sfr 
sfr_msun = np.array(data[3])
sfr = sfr_msun * 1.0e3
err_up = np.array(data[4]* 1.0e3)
err_dn = np.array(data[5]* 1.0e3)

### Make log plot
sfh_plot.bar(t_begin, sfr, width=t_width, color='w', alpha=0.6)
sfh_plot.errorbar(t_mid, sfr, yerr = [data[5]*10**3, data[4]*10**3], fmt = ' ', capsize=0, color = 'k')
sfh_plot.set_xlim(min(t_begin),10.15)
#ax.set_ylim(0,12.5)
sfh_plot.invert_xaxis()
#ax.text(0.37,0.8,'Log Age',fontsize=18)
sfh_plot.set_xlabel('Log(Time)',fontsize=22)
sfh_plot.set_ylabel('Star Formation Rate (10$^{-3}$ M$_{\odot}$ yr$^{-1}$)',fontsize=22)
sfh_plot.tick_params(axis='both', which='major', labelsize=20)

### Making the lin plot
lin_plot.plot(lin_age, csf, color='k',linestyle='solid',linewidth=3)
lin_plot.errorbar(lin_age, csf , yerr=[csf_low,csf_up], capsize=0, color='k')
lin_plot.invert_xaxis()
lin_plot.set_xlabel('Time(Gyr)',fontsize=20)
lin_plot.set_ylabel('Fraction of Stellar Mass Formed',fontsize=22)
lin_plot.tick_params(axis='both', which='major', labelsize=20)
lin_plot.set_xlim(max(lin_age), min(lin_age))
lin_plot.set_xticks(np.array([13, 12, 10, 8, 6, 4, 2, 0]))


### Making the redshift axis on top
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

ageticks = np.array([0.5, 1, 1.5, 2, 3, 6])
age_loc = 13.7 - cosmo.age(ageticks).value

rs = lin_plot.twiny()
rs.invert_xaxis()
rs.set_xlim(max(lin_age), min(lin_age))
rs.set_xticks(age_loc)
rs.set_xticklabels(['%.1f'%age for age in ageticks]);
rs.set_xlabel('Redshift (z)',fontsize = 20)
rs.xaxis.set_label_coords(0.9, 1.02)
rs.tick_params('x',labelsize=18)

### Put labels for two plots
fig.text(.14, .85, '(a)', fontsize = 22)
fig.text(.56, .85, '(b)', fontsize = 22)

plt.savefig(name +'.png')
