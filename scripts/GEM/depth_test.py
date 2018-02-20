# -*- coding: utf-8 -*-
"""
Depth test 
You should run this script in match/scripts/GEM/

Directory structure: assumes galaxy directory/sfh/ is the working directory, sfh/input_data/ should contain pars, phot and fake files
this script will create sfh/calctests/ and store results there.
Use optional -phot, -pars, -fake flags to specify their location, or change the default setting in parse_options().

Syntax: python depth_test.py [path to working directory, ending in /] -phot=[path from working directory to phot file] -pars=... -fake=...
e.g. python depth_test.py /work/04316/kmcquinn/wrangler/cusp_core/galaxies/acs/10210_UGC1281/sfh/

What it does: takes depth values from pars as a baseline
and runs calcsfh (no zinc) with 0.05 increments to +/- 0.25 from base value

Output: creates a subdirectory calctests/ and stores 120 calcsfh test runs.

After you run it, use the makeplot script (located in /work/04316/kmcquinn/wrangler/metals/scripts/plotscript/makeplot) to generate .ps files that can be viewed
"""

import subprocess as sp
from shutil import copyfile
import re
import subprocess as sp
import multiprocessing as mp
import argparse

def parse_options():
    # Creates argument parser
    parser = argparse.ArgumentParser(description='Process input parameters.')
    # Defines required arguments used on the command line
    parser.add_argument('galaxydir',action='store',help='path to galaxy directory (note: this is the whole path ending in the galaxy directory)')        
    parser.add_argument('-phot',dest='phot',action='store',help="location of phot file from galaxy directory",default='input_data/phot')
    parser.add_argument('-pars',dest='pars',action='store',help="location of pars file from galaxy directory",default='input_data/pars')
    parser.add_argument('-fake',dest='fake',action='store',help="location of fake file from galaxy directory",default='input_data/fake')
    # Parses through the arguments and saves them within the keyword args
    args = parser.parse_args()
    return args

def filtername(pars):
    '''
    opens the pars file and gets filter names
    '''
    filt=open(pars,'r')
    endl=[q[-16:-1] for q in filt if len(q)>17 and re.match('WFC...W,WFC...W',q[-16:-1])]
    red='F'+endl[0][-4:]	# red filter
    blue='F'+endl[0][3:7]	# blue filter
    filt.close()
    return blue,red

def olddepth(pars):
    '''
    pars = path to pars file, ending in .../pars
    pulls depths from pars file
    '''
    p = open(pars,'r')
    for i in range(5):
        p.readline()
    blue=p.readline().split()[:2]
    red=p.readline().split()[:2]
    depth = []
    for i in blue:
        depth.append(float(i))
    for i in red:
        depth.append(float(i))
    p.close()
    return depth
 
def makepars(pars,newpar,depth,time):
    '''
	pars = path to galaxy's pars file
	newpar = path to new parameter file
	depth = blue filter depths, red filter depths
	time = resolution being used

	uses template pars file to make a real pars file
    '''
    o = open(pars,'r')
    n = open(newpar,'w')
    for i in range(5):
        n.write(o.readline())
    bluedepth=o.readline().split()
    bluedepth[:2]=depth[:2]
    n.write(' '.join([str(i) for i in bluedepth])+"\n")
    reddepth=o.readline().split()
    reddepth[:2]=depth[2:]
    n.write(' '.join([str(i) for i in reddepth])+"\n")

    for i in range(2):
        n.write(o.readline())
    with open(time, 'r') as fobj:
        for line in fobj:
            n.write(line)

	# finished creating full pars file!
    o.close()
    n.close()

def runcalcsfh(cmdout):
    '''
    cmdout = list of calcsfh commands and output files

    this is the function to be run in parallel
    '''
    cmmd,out=cmdout
    f=open(out,"wr")
    # runs each command, waits until it's completely done, and puts the output into its corresponding output file f
    p=sp.check_call(cmmd,stdout=f)
    f.close()

def testdepth(galdir,bdep,pars,phot,fake):
    '''
    galdir = path to metals_proc directory
    bdep = depth pulled from pars file
    pars = path to galaxy's pars file
    runs depth tests!
    '''	
    testdir = galdir+'calctests/'
    sp.call(['mkdir',galdir+'calctests'])
    sp.call(['chmod','-R','770','.'])

    delta = .05             #choose how much values differ between runs
    maxdelta = .25          #choose max deviation from start values
    maxrun = int(maxdelta/delta)

    runnumber = 0 #keeps track of number of completed cycles
    runname = []
    outlist = []
    cmdlist = []
    # loops through both red and blue filters
    for i in range(-maxrun,maxrun+1):
        for j in range(-maxrun,maxrun+1):
            strrun = '%03d' % (runnumber,)
            parspath=testdir+'calcparsTEST'+strrun
            outpath=testdir+'outTEST'+strrun
            consolepath=testdir+'consoleTEST'+strrun
            tdep = bdep[:]
		 # changes depths and creates a pars file using them
            tdep[1] = bdep[1] + i * delta
            tdep[3] = bdep[3] + j * delta
            makepars(pars,parspath,tdep,'sfh_fullres')
            # creates commands to send out below and adds them to a list
            cmd=['calcsfh',parspath,phot,fake,outpath,'-Kroupa','-PARSEC']
            runname.append(runnumber)
            cmdlist.append(cmd)
            outlist.append(outpath)
            runnumber+=1
    cmmdout = []
    for i in range(len(cmdlist)):
        cmmdout.append([cmdlist[i],outlist[i]])
    # creates a pool of worker processes and runs the commands specified with map_async in parallel until it's completely done
    pool=mp.Pool()
    r=pool.map_async(runcalcsfh,cmmdout)
    pool.close()
    pool.join()


        
def main():
    # Parses through command line arguments
    args = parse_options()
    galdir = args.galaxydir
    # Defines the location of pars, phot, and fake
    pars = galdir + args.pars
    phot = galdir + args.phot
    fake = galdir + args.fake
    print filtername(pars)
    # actually runs the depth tests
    testdepth(galdir,olddepth(pars),pars,phot,fake)
    
if __name__ == "__main__":
    main()
