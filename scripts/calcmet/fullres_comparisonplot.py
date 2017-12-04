import numpy as np
from sys import argv as arg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table,Column,hstack,vstack

#Defines variables from out.final file in the table 'dat'.
def define(bop):   # bop = out.final path
    f=open(bop,'r')
    fir = f.readline()
    dats=[i.split() for i in f]
    f.close()

    # Total Star Formation Rate
    totsf = float(fir.split()[1])
    totsfu = float(fir.split()[2])
    totsfl = float(fir.split()[3])

    # Create a table with the rest of the data
    dat=Table(rows=dats,names=['start','end','dm','sfr','sfru','sfrl','met','metu','metl','uh',
                               'uhh','well','csf','csfu','csfl'],
              dtype=['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f'])
    return [[totsf,totsfu,totsfl],dat]

def getmrf(foop,nam):
        f=open(foop,'r')
        f.readline()
	pop=[i.split() for i in f]
	for k,j in enumerate(pop):
                if j[0]==nam:
                        mmrf=j[1]
                        mmrfu=j[2]
                        mmrfl=j[3]
                        break
        return [mmrf[0:4],mmrfu[0:4],mmrfl[0:4]]

def nozero(tab,time):
        ntime = []
        met = []
        metu = []
        metl = []
        for i,j in enumerate(tab['met']):
                if j != 0:
                        ntime.append(time[i])
                        met.append(j)
                        metu.append(tab['metu'][i])
                        metl.append(tab['metl'][i])
        return ntime,met,metu,metl


def main():
    zn = 'True'
    nam=arg[1]   # Galaxy directory name
    nom=arg[2]   # Galaxy name 
    cam=arg[3]   # Camera used for observation
    #zn = arg[4]  # True/False?
   
    if zn == 'True':
	zinc = ''
    elif zn == 'False':
	zinc = 'nz'

    gal='/work/04316/kmcquinn/wrangler/metals/galaxies/'+cam+'/'+nam+'/metals_proc/'

#Begin plot

    fig,axs=plt.subplots(nrows=1,ncols=2,sharex=True)

    mist=getmrf('MIST_mrflist',nom)
    pad=getmrf('PADOVA_mrflist',nom)
    par=getmrf('PARSEC_mrflist',nom)

    fig.suptitle(nom+'\n'+r'MRF$_{\mathregular{MIST}}$= $\mathregular{'+mist[0]+'^{+'+mist[1]+'}_{-'+mist[2]+'}}$, MRF$_{\mathregular{PADOVA}}$= $\mathregular{'+pad[0]+'^{+'+pad[1]+'}_{-'+pad[2]+'}}$, MRF$_{\mathregular{PARSEC}}$= $\mathregular{'+par[0]+'^{+'+par[1]+'}_{-'+par[2]+'}}$',y=1.04)
    axs[0].text(11.1,.6,'Cumulative SFH',rotation='vertical')
    axs[0].text(1.5,.6,'Age-Metallicity Relation',rotation=270)
    axs[0].text(9,-.2,'Log(Age)')
    axs[0].text(4.5,-.2,'Log(Age)')
    i = [0]
    res = ['fullres','Full']
    #for i,res in enumerate([['fullres','Full'],['no_res','N/O'],['starburst_v1res','SB V1'],['starburst_v2res','SB V2']]):

    #Define full res tables for both MIST and PADOVA, as well as saving a list of mid time between time bins
    sfm,tabm=define(gal+'sfh_fullres_MIST/out.final')
    sfpd,tabpd=define(gal+'sfh_fullres_PADOVA/out.final')
    sfpc,tabpc=define(gal+'sfh_fullres_PARSEC/out.final')
    
    time=[np.mean([tabm['start'][k],tabm['end'][k]]) for k in range(len(tabm['start']))]
    time.append(10.15)
    ax=axs[0]
    
    m=[ii for ii in tabm['csf']]
    ml=[ii for ii in tabm['csfl']]
    mu=[ii for ii in tabm['csfu']]
    m.append(0)
    ml.append(0)
    mu.append(0)
    
    pd=[ii for ii in tabpd['csf']]
    pdl=[ii for ii in tabpd['csfl']]
    pdu=[ii for ii in tabpd['csfu']]
    pd.append(0)
    pdl.append(0)
    pdu.append(0)
    
    pc=[ii for ii in tabpc['csf']]
    pcl=[ii for ii in tabpc['csfl']]
    pcu=[ii for ii in tabpc['csfu']]
    pc.append(0)
    pcl.append(0)
    pcu.append(0)


    ax.errorbar(time,m,yerr=[ml,mu],capsize=0)
    ax.errorbar(time,pd,yerr=[pdl,pdu],capsize=0)
    ax.errorbar(time,pc,yerr=[pcl,pcu],capsize=0)
    ax.set_ylim(-.1,1.1)
    ax.set_yticks([0.0,0.3,0.6,0.9])
    ax.set_xlim(10.4,6.5)
    ax.annotate(res[1], xy=(7.5, 0.05))


    mtime,mmet,mmetu,mmetl=nozero(tabm,time)
    pdtime,pdmet,pdmetu,pdmetl=nozero(tabpd,time)
    pctime,pcmet,pcmetu,pcmetl=nozero(tabpc,time)

    ax=axs[1]
    ax.errorbar(mtime,mmet,yerr=[mmetl,mmetu],fmt='.',label='MIST',capsize=0)
    ax.errorbar(pdtime,pdmet,yerr=[pdmetl,pdmetu],fmt='.',label='PADOVA ',capsize=0)
    ax.errorbar(pctime,pcmet,yerr=[pcmetl,pcmetu],fmt='.',label='PARSEC ',capsize=0)

    if max(pdmet) > max(mmet) and max(pdmet)>max(pcmet):
            big=max(pdmet)
    elif max(mmet) > max(pdmet) and max(mmet) > max(pcmet):
            big=max(mmet)
    else:
	    big=max(pcmet)
    if min(pdmet)<min(mmet) and min(pdmet) < min(pcmet):
            sma=min(pdmet)
    elif min(mmet)<min(pdmet) and min(mmet)<min(pcmet):
            sma=min(mmet)
    else:
	    sma=min(pcmet)
    
    av=abs(abs(big)-abs(sma))
    q=.3
    av2=abs(abs(big+av*q)-abs(sma-av*q))
    if av > 1:
        ax.set_yticks([round(bop,2) for bop in np.arange(sma-av2*q,big+av2*q,av2/3)])
    else:
        ax.set_yticks([round(bop,2) for bop in np.arange(sma-av2*q,big+av2*q,av2/4)])
    ax.set_ylim(sma-av*q,big+av*q)
    ax.set_xlim(10.4,6.6)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    plt.legend(handles,labels,bbox_to_anchor=(1.03,1.09),ncol=3,numpoints=1,fontsize=11)
    art=[]
    lgd=fig.subplots_adjust(hspace=0,wspace=0.1,bottom=.1,left=.1,right=.87)
    art.append(lgd)
    plt.savefig(nam+zinc+'.pdf',additional_artists=art,bbox_inches="tight")

if __name__ == '__main__':
    main()
