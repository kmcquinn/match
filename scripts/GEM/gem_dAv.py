'''
TEMPORARY INITIAL SETUP: You will have to create a directory with a name like 'sfh_fullres_MIST' (or whatever resolution/library you're using) and copy a premade parameter file into this directory until we are able to create pars files automatically with this code.

Syntax: python goodraps.py [galaxy directory name] [resolution] [library]
Optional arguments:
-phot [path to phot file from galaxy's metals_proc directory] : specifies where the phot file you would like to use is located
-fake [path to fake file from galaxy's metals_proc directory] : specifies where the fake file you would like to use is located
-nozinc : use if you do not want to include zinc in your  run
-galpath [full path to 'galaxies' directory] : specifies the main directory in which directories for all galaxy observations are located, usually just called '.../galaxies/'
'''

'''
What Goodraps Does:
- Runs initial calcsfh
- Runs hmc script to determine parameters
- pulls parameters from hmc script and runs real hmc
- runs 50 MC calcsfhs
- generates out.final files

Setup Required:
- Place phot and fake files into input_data directory located in metals_proc OR specify the paths to these files on command line with -phot and -fake arguments
- Create directory with the naming convention 'sfh_[RESOLUTION]_[LIBRARY]' in the metals_proc directory
- Create parameter file and place it in the main sfh directory (this will be automated soon!)
'''

'''
(DEVELOPMENT)
Fun things to include:
- Specify which steps to run
- Make runs with no zinc flag in separate directory
 

Does not work:
- running hmc's with no zinc flag (create nozinchmc.py which is just tony's code that works with no zinc)
- how to make new parameter file
'''

from glob import glob
import argparse
import subprocess as sp
import multiprocessing as mp
import random
import operator
import os

def parse_options():
        # Creates argument parser
        parser = argparse.ArgumentParser(description='Process input parameters.')
        # Defines required arguments used on the command line
        parser.add_argument('galaxydir',action='store',help='name of galaxy directory (note: this is just the name of the main galaxy directory, not the whole path to the directory)')
        parser.add_argument('res',action='store',help='resolution used (sfh_res/sfh_no_res/sfh_starburst_v1res/sfh_starburst_v2res')
        parser.add_argument('lib',action='store',help='stellar evolution library used for derivation (parsec/mist/padova)')
        parser.add_argument('-phot',dest='phot',action='store',help="location of phot file from galaxy's metals_proc directory",default='input_data/phot')
        parser.add_argument('-fake',dest='fake',action='store',help="location of fake file from galaxy's metals_proc directory",default='input_data/fake')
	parser.add_argument('-galpath',dest='galpath',action='store',help='full path to galaxies directory',default='/work/04316/kmcquinn/wrangler/cusp_core/galaxies/')
	parser.add_argument('-nozinc',dest='nozinc',action='store_true')
        parser.add_argument('-dAv',dest='dAv',action='store')
        # Parses through the arguments and saves them within the keyword args
        args = parser.parse_args()
        return args

# Finds the camera used for observations, flux, distance modulus, 50% completion, and magnitude depths of galaxies listed in GalCatalog (assuming it is located in same directory as oddraps)
def read_params(gal):
    params = {}
    with open('GalCatalog','r') as fobj:
			for line in fobj:
				info = line.split()
				if info[0] == gal:
					params['dir'] = info[1]
					params['flux'] = float(info[2])
					params['dist'] = float(info[3])
					params['comp'] = float(info[4])
					if len(info) > 5:
						bestDepth = info[5:9]
						params['depths'] = [float(i) for i in bestDepth]
					break
    return params

# Writes parameter file with desired inputs
def makePars(basepath, newpath, depths, times, zinc):
	#basepath: basedir
	#newpath: location, name of new pars file
	#depths: array of depth values used in this pars file
	#times: resolution the pars file should be made in; accepts 'sfh_res','sfh_no_res','sfh_starburst_v1res','sfh_starburst_v2res', or a customized time resolution file in oddraps directory
	
	#open base pars file and write first 5 lines to new pars file
	g = open(basepath+'input_data/pars', 'r')
	p = open(newpath, 'w')
	p.write(g.readline())
	line2 = g.readline()[:-1]
	if zinc == True:
 ########### GET REAL ZINC VALUES HERE ##############
		line2 = ZINCFILES
	else:
		line2 = "-2.0 0.1 0.1"
	line2 = line2+"\n"
	p.write(line2)
	for i in range(0,3):
		p.write(g.readline())
	
	#now we are at the filter lines
	blueline = g.readline().split()
	blueline[:2] = depths[:2]	#replace depths found in pars file with desired
	p.write(' '.join([str(n) for n in blueline])+"\n")
	
	redline = g.readline().split()
	redline[:2] = depths[2:]
	p.write(' '.join([str(n) for n in redline])+"\n")
	
	#next two lines can be copied over
	p.write(g.readline())
	p.write(g.readline())
	
	#now we need to add the requested timebin
	with open(times, 'r') as fobj:
		for line in fobj:
			p.write(line)
			
	#desired pars file has now been created
	p.close()
	g.close()

# runs initial calcsfh
def run_initial_calcsfh(mp_path,sfh_path, args):
	photLoc = mp_path+args.phot
	fakeLoc = mp_path+args.fake
        dAv = args.dAv
	#run calcsfh once for use with hybridMC
	######## CHANGE PATH TO PARAMETER FILE ONCE CODE CAN GENERATE PARAMETER FILES ###################
	comm1 = "calcsfh "+sfh_path+"pars "+photLoc+" "+fakeLoc+" "+sfh_path+"out -Kroupa "+"-dAv="+args.dAv+" "
	if (args.lib == "MIST") or (args.lib == "PARSEC"):
		comm2 = "-"+args.lib.upper()+" "
	else:
		comm2 = " "
	comm3 = "-mcdata"
	if args.nozinc == True:
		comm4 = ""
	else:
		comm4 = " -zinc"
	comm = comm1 + comm2 + comm3 + comm4
	f = open(sfh_path+"console.txt", "wb")
	sp.call(comm.split(),stdout=f)
	f.close()

# runs HMC
def run_hmc(sfh_path):
	#run code to find hmc parameters
	comm = "python hmc.py "+sfh_path
	f = open(sfh_path+'hmc_console_output.dat','w')
	sp.call(comm.split(),stdout=f)
	f.close()
	#pull hmc parameters from output
	f = open(sfh_path+'hmc_param_console.dat')
	pars = [i.split() for i in f]
	f.close()
	# hybridMC command to be run should be last line in hybrid_params.dat file
	comm = pars[-1]
	#run hybridMC 1000 times
	f = open(sfh_path+"hybrid_console.txt", "w")
	sp.call(comm,stdout=f)
	f.close()

# dummy function for parallelizing 50 MC calcsfh's
def Calcwork(arr):
	comm, output = arr[0], arr[1]
	f = open(output, "wb")	
	sp.call(comm.split(),stdout=f)
	f.close()
	return 0

# runs 50 MC calcsfh's
def run_calcmc(mp_path,sfh_path,args,params):
	photLoc = mp_path+args.phot
	fakeLoc = mp_path+args.fake
	goodDepths = params['depths']
	#need to figure out values for logterrsig, mbolerrsig
	g = open(sfh_path+"pars", 'r')
	Mt = float(g.readline().split()[0])
	g.close()
	Mv = goodDepths[1]
	Mi = goodDepths[3]
	Mv_diff = Mv - Mt
	Mi_diff = Mi - Mt
	#pick choice that minimizes total error in *_diff
	#teffdata in Desktop dir has it all
	err = []
	vals = []
	with open('teffdata', 'r') as fobj:
		fobj.readline()
		for line in fobj:
			row = line.split()
			vals.append([float(row[2]), float(row[3])])
			err.append(abs(Mv_diff - float(row[0])) + abs(Mi_diff - float(row[1])))
	minloc, minval = min(enumerate(err), key=operator.itemgetter(1))
	lgsig, mbol = vals[minloc][0], vals[minloc][1]
	part1 = "calcsfh "+sfh_path+"pars "+photLoc+" "+fakeLoc+" "+sfh_path+"out_"
	part2 = " -Kroupa "
	if (args.lib == "MIST") or (args.lib == "PARSEC"):
		part3 = "-"+args.lib+" "
	else:
		part3 = ""
	if args.nozinc == True:
		part4 = "-mcdata -mcseed="
	else:
		part4 = "-zinc -mcdata -mcseed="
		
	part5 = " -logterrsig="
	part6 = " -mbolerrsig="
	
	runarr = []
	for i in range(0,50):
		rand = random.randint(0,5000)
		digit = '%02d' % (i,)
		outfile = open(sfh_path+"mcseed_values.txt","a")#Storing the values of mcseed.
		out_string = "mcseed value " + str(digit) + "=" + str(rand) + "\n"
		outfile.write(out_string)
		outfile.close()
		comm = part1 + digit + part2 + part3 + part4 + str(rand) + part5 + str(lgsig) + part6 + str(mbol)
		runarr.append([comm,sfh_path+"console"+digit])
	pool = mp.Pool(None)
	pool.map_async(Calcwork, runarr)	#fills cpu cores with dowork jobs, each with different flail value from runname
	pool.close()
	pool.join()

# Takes path to directory that contains all of the previous runs (initial calcsfh, hmc, 50 MC calcsfh's) and combines the results to generate a final results file.
def generate_final(fullpath,galname):      # fullpath = sfh_path
    comarr = []
    comarr.append("zcombine "+fullpath+"out > "+fullpath+"out.zc")
    comarr.append("zcombine -unweighted -medbest -jeffreys -best="+fullpath+"out.zc "+fullpath+"out.mcmc > "+fullpath+"out.mcmc.zc")
    comarr.append("zcombine -unweighted -medbest -best="+fullpath+"out.zc "+fullpath+"out_?? > "+fullpath+"out.sys.zc")
    comarr.append("zcmerge "+fullpath+"out.zc "+fullpath+"out.mcmc.zc "+fullpath+"out.sys.zc -absolute > "+fullpath+"out.final")
    #combining wih hybrid uncertainities for comparison
    comarr.append("zcmerge "+fullpath+"out.zc "+fullpath+"out.mcmc.zc -absolute > "+fullpath+"out.hybrid.final")
    comarr.append('python 9_panel.py '+fullpath+'out.cmd '+fullpath+'out.final '+fullpath+galname+'.png')
    for i in comarr:
	os.system(i)
	#so we have created the pars file used in the main calcsfh runs, and completed a  calcsfh analysis of this galaxy
	#a plot has been created showing the fit and uncertainties

def main():
    # Parses through command line arguments
    args = parse_options()
    # Extracts data for specified galaxy from GalCatalog
    params = read_params(args.galaxydir)

    # Ensures the input path to the full 'galaxies' directory ends in a slash to avoid future error
    if args.galpath[-1] != '/':
	args.galpath = args.galpath + '/'

    # Specifies path to metals_proc directory within main galaxy directory
    mp_path = args.galpath+params['dir']+"/"+args.galaxydir+"/sfh/"
    inpath = mp_path + 'input_data/'
    sfh_path = mp_path+args.res+'_'+args.lib.upper()+'/'

    # if the sfh_[res]_[lib] directory does not exist in basepath, creates it
    if args.res+'_'+args.lib.upper() not in glob(mp_path+'*'):
    	comm = "mkdir "+sfh_path
#    	sp.call(comm.split())

    run_initial_calcsfh(mp_path,sfh_path,args)
    run_hmc(sfh_path)
    run_calcmc(mp_path,sfh_path,args,params)
    generate_final(sfh_path,args.galaxydir)


if __name__ == "__main__":
    main()
