import subprocess as sp
import sys
from shutil import copyfile
import re
import subprocess as sp
import multiprocessing as mp

## HOW TO RUN DTEST ##
# python dtest.py [name of galaxy directory] -phot=[path to phot file from galaxy's metals_proc directory] -pars=[path to parameter file from metals_proc] -fake=[path to fake file from metals_proc] -run=[all/fit/scriptdir]
# phot, pars, and fake arguments are all optional arguments that help set up the scriptdir directory
# run is also an optional argument that allows certain parts of the script to be run
  # accepted arguments:
    # all (default): sets up scriptdir and runs through depth tests
    # fit : runs through depth tests
    # setup : sets up scriptdir


def setup(bdir,par=None,pho=None,fa=None):
	# bdir = galaxy's metals_proc directory
	# par = path to pars file from galaxy's metals_proc directory
	# pho = path to phot file from galaxy's metals_proc directory
	# fa = path to fa file from galaxy's metals_proc directory
	
	# makes scriptdir directory, if it doesn't exist already
	if "scriptdir" in sp.check_output(["ls",bdir]).splitlines():
		print "scriptdir exists"
	else:
		sp.call(["mkdir",bdir+"scriptdir/"])

	# find pars file in conf_new_dol if none specified
	if par==None:
		pardir=bdir+"../conf*/"
		parfile=sp.check_output(["ls",pardir]).splitlines()
		parfile=[i for i in parfile if 'matchpars' in i]
		if len(parfile)==1: # if there is only one parameter file in conf_new_dol, uses this one
			copyfile(pardir+parfile[0],bdir+"scriptdir/basepars")
			print "found single parameter file in conf* directory"
		else: # for more than one parameter file
			parpar= [i for i in parfile if i[-12:]=='matchpars231']
                        if len(parpar)==1: # checks for *matchpars231 and if only one file is found, uses this
				copyfile(pardir+parpar[0],bdir+"scriptdir/basepars")
			else: # if multiple *matchpars231 files, stops setFolder
				print "ERROR: Multiple parameter files found; please specify which to use"
				return None
	# if pars file is spedified, copies it to scriptdir
	else:
		copyfile(bdir+par,bdir+"scriptdir/basepars")
	
	# pulls filter values from parameter file
	pf=bdir+'scriptdir/basepars'
	filt=open(pf,'r')
	endl=[q[-16:-1] for q in filt if len(q)>17 and re.match('WFC...W,WFC...W',q[-16:-1])]
	red='F'+endl[0][-4:]	# red filter
	blue='F'+endl[0][3:7]	# blue filter
	filt.close()
	
	# copies phot file from proc_new_dol if none spedified
	if pho==None:
		phodir=bdir+'../proc*/'
		phofile=sp.check_output(["ls",phodir]).splitlines()
		# looks for phot file with same filters specified in pars file
		phofile=[i for i in phofile if '.gst.match' in i and red in i and blue in i]
		copyfile(phodir+phofile[0],bdir+'scriptdir/phot')
	
	# if phot file specified, copies it to scriptdir
	else:
		copyfile(bdir+pho,bdir+'scriptdir/phot')	
	
	# finds and copies fake file in proc_new_dol if none specified
	if fa==None:
		fadir=bdir+'../proc*/'
		fafile=sp.check_output(["ls",phodir]).splitlines()
                fafile=[i for i in fafile if 'gst.matchfake' in i]
		copyfile(fadir+fafile[0],bdir+'scriptdir/fake')
	# if fake file specified, copies it to scriptdir
	else:
		copyfile(bdir+fa,bdir+'scriptdir/fake')
	# returns the filters in the form [F___,F___]
	return blue,red


def basepars(sdir):	# creates template pars file
	# sdir = path to scriptdir

	o = open(sdir+"basepars",'r')
	old= [i for i in o]
	o.close()
	o = open(sdir+"basepars",'r')
	n = open(sdir+'pars','w')
	header = o.readline().split()   #contains IMF (possibly), dist, extinction

        # get rid of IMF
        if len(header) == 7:
        	nheader = header[1:]    #new header skips first entry
        else:
        	nheader = header
	# write out new header in same row format
	n.write(' '.join([str(i) for i in nheader])+"\n")
	# write correct zinc line
	n.write('-2.0 0.1 0.1\n')
	# copies line from old file
	n.write(old[2])
	# red clump line; just needs the number 1 for now
	n.write('1\n')
	# copies the next 3 filter lines
	n.write(old[4])
	n.write(old[5])
	n.write(old[6])
	# blank line, and two zeroes underneath
	n.write('\n')
	n.write('0 0\n')
	# finished creating pars file!
	n.close()
	o.close()

def makepars(bdir,newpar,depth,time):
	# bdir = path to galaxy's metals_proc directory
	# newpar = path to new parameter file
	# depth = blue filter depths, red filter depths
	# time = resolution being used

	# uses template pars file to make a real pars file
	o = open(bdir+'scriptdir/pars','r')
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

def olddepth(pdir):
	# pdir = path to basepars file

	# pulls depths from basepars file
	p = open(pdir,'r')
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


def runcalcsfh(cmdout):
	# cmdout = list of calcsfh commands and output files

	# this is the function to be run in parallel
	cmmd,out=cmdout
	f=open(out,"wr")
	# runs each command, waits until it's completely done, and puts the output into its corresponding output file f
	p=sp.check_call(cmmd,stdout=f)
	f.close()
	# sp.check_call(['pg_cmd',out+'.cmd',out+'.ps'])

def testdepth(bdir,bdep):
	# bdir = path to metals_proc directory
	# bdep = depth pulled from pars file
	
	# runs depth tests!
	scrdir=bdir+'scriptdir/'
	testdir=scrdir+'calctests/'
	sp.call(['mkdir',scrdir+'calctests'])
	phot=scrdir+'phot'
	fake=scrdir+'fake'

        delta = .05             #choose how much values differ between runs
        maxdelta = .25          #choose max deviation from start values
        runnum = 1              #keeps track of number of completed cycles
        maxrun = int(maxdelta/delta)

	runnumber = 0
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
			makepars(bdir,parspath,tdep,'sfh_fullres')
			# creates commands to send out below and adds them to a list
			cmd=['calcsfh',parspath,phot,fake,outpath,'-Kroupa','-PARSEC']
			runname.append(runnumber)
			cmdlist.append(cmd)
			outlist.append(outpath)
			runnumber+=1
	cmmdout = []
	for i in range(len(cmdlist)):
		cmmdout.append([cmdlist[i],outlist[i]])
	# creates a pool of worker processes and runs the commands specified wiht map_async in parallel until it's completely done
	pool=mp.Pool()
	r=pool.map_async(runcalcsfh,cmmdout)
	pool.close()
	pool.join()


def main():
	GalName= sys.argv[1]
        with open('GalCatalog','r') as fobj:
                        for line in fobj:
                                info = line.split()
                                if info[0] == GalName:
                                        GalDir = info[1]
					break
	basedir = "/work/04316/kmcquinn/wrangler/metals/galaxies/"+GalDir+"/"+GalName+"/metals_proc/"
	scriptdir = basedir+"scriptdir/"
	
	phot=None
	fake=None
	pars=None
	run='all'

	for i in sys.argv:
		if i[:4]=='phot':
			phot=i[5:]
		if i[:4]=='fake':
			fake=i[5:]
		if i[:4]=='pars':
			pars=i[5:]
		if i[:3]=='run':
			run=i[4:]

	if run=='all':
		print setup(basedir,par=pars,pho=phot,fa=fake)
		basepars(scriptdir)
		testdepth(basedir,olddepth(scriptdir+'basepars'))
	elif run == 'setup':
                print setup(basedir,par=pars,pho=phot,fa=fake)
                basepars(scriptdir)
                testdepth(basedir,olddepth(scriptdir+'basepars'))
	elif run=='fit':
		testdepth(basedir,olddepth(scriptdir+'basepars'))
	else:
		print "run = "+run

if __name__ == "__main__":
	main()
