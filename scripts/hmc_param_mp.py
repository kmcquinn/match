#Written by Anthony Pahl 6/12/17
import sys
import numpy as np
import multiprocessing as mp
import subprocess32 as subprocess
import time
import shlex
import os
import glob

##hmc_param_mp.py
##Finds dt, tint, and sscale parameters for NON-ZINC runs of hybridMC.
##Uses parallel programming to average short hMC runs to quickly
##generate these parameters by the process described by the README of 
##match v2.7. 

##Input: output files from calcsfh
##mine are in the format cout_*.dat
datfs=glob.glob('cout_*.dat')

##Output: a console file that summarizes hmc_param_mp's activity
consolefile = 'hmc_param_console.dat'
if os.path.isfile(consolefile):
    os.remove(consolefile)
##and an output file that contains properly formatted hMC runs with 
##the correct parameters
cmdsfile = 'hmc_param_out.dat'
if os.path.isfile(cmdsfile):
    os.remove(cmdsfile)

##Set Number of hMC runs to average when determining % acceptance or 
##ratio of sigmas from first step to last step.
nhmcs = 8




def hMC_wrapper(args):
    result = hMC(*args)
    #check for timeouts
    try:
        return result.communicate(timeout=300)[0]
    except subprocess.TimeoutExpired:
        result.terminate()
        return np.nan
        
def hMC(ifile,ofile,flags,seed):
    cmd= ' '.join([hybridmc,ifile,ofile,flags,'-mcseed='+str(seed)])
    return subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE)

def get_accept_frac(sresult):
    if np.array(sresult).dtype == float and np.isnan(sresult):
        return np.NaN
    stemp = sresult.split('\n')
    if (stemp[-1]):
        raise Exception('result string: something went wrong')
    stemp = stemp[-2]
    frac = stemp.split()[1].split('/')
    return float(frac[0]) / float(frac[1])

def get_sig_frac(sresult):
    if np.array(sresult).dtype == float and np.isnan(sresult):
        return np.NaN
    stemp = sresult.split('\n')
    if (stemp[-1]):
        raise Exception('result string: something went wrong')
    ftemp = stemp[-2]
    st=ftemp.find('SD = ')+5
    end=ftemp.rfind(')')
    final_sig = float(ftemp[st:end])

    itemp = stemp[2]
    st=itemp.find('SD = ')+5
    end=itemp.rfind(')')
    init_sig = float(itemp[st:end])

    return init_sig/final_sig
    
     
try:
    hybridmc = subprocess.check_output(['which','hybridMC']).rstrip()
except subprocess.CalledProcessError as err:
    raise Exception("hybridMC not found, try loading match module.")

mcmcfs=[]
mcmcfs_out=[]
hconsolefs=[]
for i in datfs:
    nstr = i.lstrip('cout_')
    nstr = nstr.rstrip('.dat')
    mcmcfs.append('hmc_param_'+nstr+'.mcmc')
    mcmcfs_out.append('cout_'+nstr+'.mcmc')
    hconsolefs.append('hybrid_console_'+nstr+'.dat')



cpus = mp.cpu_count()
p = mp.Pool(processes = cpus)

n_async = cpus / nhmcs #number of steps to run asynchr.
if n_async == 0:
    n_async = 1
out_sh = '' #string with output hybridMC commands

with open(consolefile,'a') as consolef:
    print >>consolef,'CPUS: ',cpus
    print >>consolef,'HMCS TO AVERAGE: ',nhmcs
    print >>consolef,'HYBRIDMC TO USE: ',hybridmc

#for number of inputs...
for i,e in enumerate(datfs):
    inputf = datfs[i]
    outputf = mcmcfs[i]
    outputf_cmd = mcmcfs_out[i]
    hconsolef = hconsolefs[i]
    with open(consolefile,'a') as consolef:
        print >>consolef,'INPUT: ',inputf
    #initial guesses
    nmc_ = 10
    nsfr_ = 0
    dt_ = 0.01
    tint_ = 0.01
    a_frac = 1.

    #determine dt
    while (a_frac > 0.9):
        results_arr = np.empty(n_async,dtype=object)
        dt_t = dt_
        tint_t = tint_
        for j in range(n_async):
            dt_t+=0.005
            tint_t+=0.005
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_t),'-tint='+str(tint_t)])
            seeds = np.random.randint(0,10000,nhmcs)
            iter_args = [(inputf,outputf,flags,)]*nhmcs
            iter_args = [e+(seeds[i],) for i,e in enumerate(iter_args)]
            results_arr[j] = p.map_async(hMC_wrapper,iter_args)
        for map_results in results_arr:  
            results = map_results.get()
            dt_+=0.005
            tint_+=0.005
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),'-tint='+str(tint_)])
            a_frac = np.average([get_accept_frac(i) for i in results])
            with open(consolefile,'a') as consolef:
                print >>consolef,flags,' => ',a_frac,' accepted'
            if (a_frac < 0.9):
                break
            
    dt_-=0.005
    tint_-=0.005
    
    #determine tint
    sig_frac = 0.
    while (sig_frac < (2./3.)):
        results_arr = np.empty((n_async),dtype=object)
        tint_t = tint_
        for j in range(n_async):
            tint_t+=dt_
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),'-tint='+str(tint_t)])
            seeds = np.random.randint(0,10000,nhmcs)
            iter_args = [(inputf,outputf,flags,)]*nhmcs
            iter_args = [e+(seeds[i],) for i,e in enumerate(iter_args)]
            results_arr[j] = p.map_async(hMC_wrapper,iter_args)
        for map_results in results_arr:  
            results = map_results.get()
            tint_+=dt_
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),'-tint='+str(tint_)])
            sigs = [get_sig_frac(i) for i in results]
            sig_frac = np.average(sigs)
            if 0 in sigs:
                zeros = np.where(np.array(sigs) == 0.)
                sigs_no_z = np.delete(sigs,zeros)
                if sigs_no_z.size ==0:
                    av_noz = 0.
                else:
                    av_noz = np.average(sigs_no_z)
                zerostr = ', '+str(len(zeros[0]))+'/'+str(len(sigs))+' zeros found (with: '+str(sig_frac)+')'
                sig_frac=av_noz
            else:
                zerostr = ''
            with open(consolefile,'a') as consolef:
                print >>consolef,flags,' => ',sig_frac,' ratio of sigmas',zerostr
            #end = time.time() ##
            #print end-start,' s' ##
            if (sig_frac > (2./3.)):
                break

    #determine sscale
    nsfr_=10
    sscale_=0.
    a_frac = 0.
    while (a_frac < 0.80 or np.isnan(a_frac)):
        results_arr = np.empty(n_async,dtype=object)
        sscale_t = sscale_
        for j in range(n_async):
            sscale_t+=1.
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),'-tint='+str(tint_),
                              '-sscale='+str(sscale_t)])
            seeds = np.random.randint(0,10000,nhmcs)
            iter_args = [(inputf,outputf,flags,)]*nhmcs
            iter_args = [e+(seeds[i],) for i,e in enumerate(iter_args)]
            results_arr[j] = p.map_async(hMC_wrapper,iter_args)
        for map_results in results_arr:  
            sscale_+=1.
            results = map_results.get()
            fracs=[get_accept_frac(i) for i in results]
            a_frac = np.average(fracs)
            if np.isnan(np.min(fracs)):
                nans = [np.isnan(i) for i in fracs]
                fracs_noz = np.array([i for i in fracs if not np.isnan(i)])
                if fracs_noz.size == 0:
                    av_noz = np.nan
                else:
                    av_noz = np.average(fracs_noz)
                nanstr = ', '+str(len(nans))+'/'+str(len(fracs))+' timeouts occured (without: '+str(av_noz)+')'
            else:
                nanstr = ''
            flags = ' '.join(['-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),'-tint='+str(tint_),
                              '-sscale='+str(sscale_)])
            with open(consolefile,'a') as consolef:
                print >>consolef,flags,' => ',a_frac,' accepted'+nanstr
            if not (a_frac < 0.80 or np.isnan(a_frac)):
                break

    nmc_=100
    nsfr_=100
    sburn_=10

    out_sh_t = 'hybridMC '+' '.join([inputf,outputf_cmd,'-nmc='+str(nmc_),'-nsfr='+str(nsfr_),'-dt='+str(dt_),
                                   '-tint='+str(tint_),'-sscale='+str(sscale_),'-sburn='+str(sburn_)])+ ' > '+hconsolef
    out_sh += out_sh_t + '\n'
    with open(consolefile,'a') as consolef:
        print >>consolef,out_sh_t
    with open(cmdsfile,'a') as cmdsf:
        print >>cmdsf,out_sh_t

    #someone clean up these damn files
    os.remove(outputf)
