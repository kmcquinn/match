import time
import subprocess as sp
import sys
import random
import math
import operator
import commands
from shutil import copyfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fnmatch
'''
This version is just going to get us off the ground in running analysis
It sends every command out linerally. This includes depth optimization and hybridMC runs, which would stand the greatest improvement from parallel processing.
oddraps = on dwarf disks, running a python script

Ensure script has all files needed in same folder:
	GalCatalog
	sfh_* files
	batch script to run on TACC
	
Edit batch file FitBatch to run this script (command: python FitTesting.py)
	set time accordingly (currently about 30 hours for full run)
	run through login node: "sbatch FitBatch"

Failure will most probably occur in the setFolder function
	Python error output will be shown in the console output set in FitBatch
	Assumes the given calcsfh files and pars file are in the folders shown
	Can probably edit this to point to the right places
'''

def setFolder(bir):
#passed
	'''
	hey, this guy is taking the first file it sees that fits those wildcards
	that means if there is more than one match, any of those files is fair game for it to pull
	sorry
	'''
	#this creates a directory in the work folder and moves the given files (pars,phot,fake) over for safety
	#first setup folder in gal dir to build files
	sp.call(["mkdir",bir+"scriptdir/"])
	#find pars file in input_data folder
	fold = bir+"../conf_new_dol/"
	inList = sp.check_output(["ls",fold]).splitlines()
	cull = fnmatch.filter(inList,"*matchpars231")
        copyfile(fold+cull[0], bir+"scriptdir/basepars")
	#find phot file in input_data folder
	#NOTE: Here I will also store the written filter values for editFiles to use
	fold = bir+"../proc_new_dol/"
	inList = sp.check_output(["ls",fold]).splitlines()
	exten = "*.gst.match"
	cull = fnmatch.filter(inList,exten)
	#store filter values
	FiltA = cull[0][len(cull)-len(exten)-11:len(cull)-len(exten)-6]
	FiltB = cull[0][len(cull)-len(exten)-5:len(cull)-len(exten)]
	#copy phot file to script dir
	copyfile(fold+cull[0], bir+"scriptdir/phot")
	#find fake file in input_data folder
	cull = fnmatch.filter(inList,"*gst.matchfake")
	#copy fake file to script dir
	copyfile(fold+cull[0], bir+"scriptdir/fake")
	return FiltA, FiltB
	#ok so all the files we need are in our new directory, and we grabbed the filter names from the filenames to compare with the given pars file. 

def editFiles(sir, afile, bfile):
	#now we need to make any needed changes to the pars file before using them in calcsfh
	#phot: good to go
	#fake: good to go
	#pars: could include IMF, false filter info, wrong timebin, etc.
	g = open(sir+"basepars", 'r')	#open given pars
	p = open(sir+"pars", 'w')		#create edited calcsfh pars
	header = g.readline().split()	#contains IMF (possibly), dist, extinction
	
	#get rid of IMF
	if len(header) == 7:
		nheader = header[1:]	#new header skips first entry
	else:
		nheader = header
	
	p.write(' '.join([str(n) for n in nheader])+"\n")	#write out new header in same row format
	#now change second row if needed
	secline = g.readline().split()
	if float(secline[0]) != -2.2:
		secline[0] = -2.2
	p.write(' '.join([str(n) for n in secline])+"\n")
	
	#don't care about these next two rows
	p.write(g.readline())
	p.write(g.readline())
	
	#make sure filter names are correct
	filine = g.readline().split()
	filt = filine[-1]
	comm = filt.find(',')
	filtA = filt[:comm]
	filtB = filt[comm+1:]
	if filtA[3:] != afile[1:] or filtB[3:] != bfile[1:]:	#check to see if pars filters match file filters
		#file name is going to have superiority here
		filtA[3:] = afile[1:]
		filtB[3:] = bfile[1:]
	filtstring = str(filtA)+','+str(filtB)
	filine[-1] = filtstring
	p.write(' '.join([str(n) for n in filine])+"\n")
	
	#next two lines for filter depths and names
	#read line and replace filter name with one used in previous line
	blueline = g.readline().split()
	blueline[-1] = filtA
	startblue = blueline[:2]	#contains starting blue depths
	p.write(' '.join([str(n) for n in blueline])+"\n")
	redline = g.readline().split()
	redline[-1] = filtB
	startred = redline[:2]		#contains starting red depths
	p.write(' '.join([str(n) for n in redline])+"\n")
	
	#next line should be empty, must check to make sure
	blankline = g.readline()
	if len(blankline) != 1:
		p.write('\n')	#in this case blankline has useful data, need to write it
		p.write(blankline)
	else:
		p.write("\n")
		p.write(g.readline())	#must read another line to catch up with the first case
	
	#here we have the number of timebins and timebins left in pars
	#the plan is to leave this pars file as is, and construct the rest of it with the timebin files based on the calcsfh run it's used for
	p.close()
	g.close()
	return [float(i) for i in startblue] + [float(i) for i in startred]	#returns startblue and startred as one list of float entries

def makePars(basepath, newpath, depths, times):
	#produced calcsfh pars file based on given parameters
	#basepath: basedir
	#newpath: location, name of new pars file
	#depths: array of depth values used in this pars file
	#times: which timebin to append pars file with, assume this is a string with the name of the timebin file in /galaxies/ dir
	
	#open base pars file and write first 5 lines to new pars file
	g = open(basepath+'scriptdir/pars', 'r')
	p = open(newpath, 'w')
	for i in range(0,5):
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



def runBatch(comm, name, times):
	#creates batch script w/ given command and runs, waiting for completion
	#first need to make the file
	p = open("basebatch.txt", "r")
	g = open("script","w")
	g.write(p.readline())
	g.write(p.readline()[:-1] + name + "\n")
	g.write(p.readline()[:-1] + name + ".o%j" + "\n")
	g.write(p.readline()[:-1] + "1" + "\n")
	g.write(p.readline())
	g.write(p.readline()[:-1] + times + "\n")
	g.write(p.readline())
	g.write(comm)
	p.close()
	g.close()
	
	#now need to queue script and find job id used to check status
	jobnum = int(sp.check_output(["sbatch","script"]).splitlines()[-1].split()[-1])
	
	#check for job completion every 30 seconds
	'''
	here i'm going to have to be sloppy, as it seems like squeue is pretty unreliable
	at showing completed jobs
	i'm just going to assume a job is complete if it doesn't say pending or running
	i'm sure this will work well
	'''
	while True:
		line = sp.check_output(["squeue","-j",str(jobnum)]).splitlines()[-1].split()
		if (line[4] != "R") and (line[4] != "PD"):
			break
		time.sleep(30)


def grabDepths(outpath):
	#only for calcsfh pars files
	deparr = []
	g = open(outpath, 'r')
	for i in range(0,5):
		g.readline()
	bluefilt = g.readline().split()[:2]
	for x_i in bluefilt:
		deparr.append(float(x_i))
	redfilt = g.readline().split()[:2]
	for x_i in redfilt:
		deparr.append(float(x_i))
	g.close()
	return deparr
	#returns four filter depth values


def useLauncher(name,commarr, time):
	#makes launcher script with given commands, runs and waits for completion
	#first make launcher.slurm file
	g = open("pylaunch.slurm","w")
	#header info
	with open("launchhead.txt","r") as fobj:
		for line in fobj:
			g.write(line)
	#body info
	b = open("launchbody.txt","r")
	g.write(b.readline()[:-1] +name+"\n")
	g.write(b.readline()[:-1] +str(math.ceil(float(len(commarr))/16.0))+"\n") #one core per command
	g.write(b.readline()[:-1] +str(len(commarr))+"\n")
	g.write(b.readline())
	g.write(b.readline()[:-1] +name+".o%j"+"\n")
	g.write(b.readline()[:-1] + time + "\n")
	b.close()
	#footer info
	with open("launchfoot.txt","r") as fobj:
		for line in fobj:
			g.write(line)
	g.close()
	
	#now need to build paramlist
	p = open("paramlist","w")
	for i in comarr:
		if i == len(comarr) - 1:
			p.write(i)
			continue
		p.write(i + "\n")
	p.close()
	
	#now launch slurm script, grab job id for tracking
	jobnum = int(sp.check_output(["sbatch","pylaunch.slurm"]).splitlines()[-1].split()[-1])
	
	#check for job completion every 30 seconds
	'''
	here i'm going to have to be sloppy, as it seems like squeue is pretty unreliable
	at showing completed jobs
	i'm just going to assume a job is complete if it doesn't say pending or running
	i'm sure this will work well
	'''
	while True:
		line = sp.check_output(["squeue","-j",str(jobnum)]).splitlines()[-1].split()
		if (line[4] != "R") and (line[4] != "PD"):
			break
		time.sleep(30)

def calcFit(bir, scr, filtStart):
	'''
	this is a modified version of the calcsfh depth optimatization on oddraps
	it will take the starting filter values in the given par file and try every variation of the weak V and I filter limits witin a given range and step size
	Each CMD's filter and resulting fit is given in the FiltResults file located in the calctests folder
	'''
	#runs simple calcsfh many times to find filter depths that produce best fit
	#create folder to house calcsfh test runs
	sp.call(["mkdir",scr+"calctests"])
	scrstring = scr+"calctests/"
	scrPhot = scr+"phot"	#save path to phot and fake files for convienence
	scrFake = scr+"fake"
	
	#For now, not going to incorporate random depth runs. It would probably go here anyway
	#Need to grab initial filter values from the editFiles function
	sold = []	#stores base fit, start depths
	sold.append(0)
	sold = sold + filtStart
	
	#before starting all the runs, I need to run calcsfh once with these starting depths to get a baseline to compare all the permutations with
	#first step is to create the pars file with the starting depths
	makePars(bir, scrstring+"parsFirst", sold[1:], "sfh_fullres")
	comm = "calcsfh "+scrstring+"parsFirst "+scrPhot+" "+scrFake+" "+scrstring+"outFirst -Kroupa -PARSEC"
	fitstr = sp.check_output(comm.split()).splitlines()[-1].split()[-1]
	#check consolefile to record starting fit value
	#format: 'fit=3573.515152'
	fitval = float(fitstr[4:])
	print("fitval is: "+str(fitval))
	#record this value to sold array
	sold[0] = fitval
	delta = .05		#choose how much values differ between runs
	maxdelta = .25		#choose max deviation from start values
	runnum = 1		#keeps track of number of completed cycles
	maxrun = int(maxdelta/delta)
	
	flail = 0
	fitlist = []
	permlist = []
	g = open(scrstring+"FiltResults","w")
	g.write("Run Number\tDepth1\tDepth2\tFit Value\n")
	for i in range(-maxrun,maxrun+1):
		for j in range(-maxrun,maxrun+1):
			strflail = '%03d' % (flail,)		#convert run number to string for out files
			parspath = scrstring+"calcparsTEST"+strflail
			outpath = scrstring+"outTEST"+strflail
			consolepath = scrstring+"consoleTEST"+strflail
			totest = sold[1:]
			totest[1] = totest[1] + i * delta
			totest[3] = totest[3] + j * delta
			makePars(bir, parspath, totest, "sfh_fullres")	#make pars file with 'test' to indicate temp file
			comm = "calcsfh "+parspath+" "+scrPhot+" "+scrFake+" "+outpath+" -Kroupa -PARSEC"
			fitout = float(sp.check_output(comm.split()).splitlines()[-1].split()[-1][4:])
			sp.call(["pg_cmd",outpath+".cmd",outpath+".ps"])
			fitlist.append(fitout)
			permlist.append(totest)
			g.write(strflail+"\t"+str(totest[1])+"\t"+str(totest[3])+"\t"+str(fitout)+"\n")
			flail = flail + 1
	minloc, minval = min(enumerate(fitlist), key=operator.itemgetter(1))
	beststr = '%03d' % (minloc,)
	g.write("Best run: outTEST"+beststr+", "+str(minval)+" with filter values "+str(permlist[minloc][1])+" "+str(permlist[minloc][3]))
	g.close()
	return permlist[minloc]

def fullCalc(bpath, fullpath, goodDepths, tbins):
	#runs the full calcsfh workflow, incudes hybridMC, .ps plot of results
	photLoc = bpath+"scriptdir/phot"
	fakeLoc = bpath+"scriptdir/fake"
	#create needed pars file, run full calcsfh script
	makePars(bpath, fullpath+"pars", goodDepths, tbins)
	#need to figure out values for logterrsig, mbolerrsig
	g = open(fullpath+"pars", 'r')
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
##This is a bad way to calculate error
	with open('teffdata', 'r') as fobj:
		for line in fobj:
			row = line.split()
			vals.append([float(row[2]), float(row[3])])
			err.append(abs(Mv_diff - float(row[0])) + abs(Mi_diff - float(row[1])))
	minloc, minval = min(enumerate(err), key=operator.itemgetter(1))
	lgsig, mbol = vals[minloc][0], vals[minloc][1]
	#run calcsfh once for use with hybridMC
	comm = "calcsfh "+fullpath+"pars "+photLoc+" "+fakeLoc+" out -Kroupa -PARSEC -mcdata > "+fullpath+"console.txt"
	sp.call(comm.split())
	#run hybridMC 1000 times
	comm = "hybridMC "+fullpath+"out.dat "+fullpath+"out.mcmc -tint=2.0 -nmc=10000 -dt=0.015 > "+fullpath+"hybrid_console.txt"
	sp.call(comm.split())
	
	part1 = "calcsfh "+fullpath+"pars "+photLoc+" "+fakeLoc+" "+fullpath+"out_"
	part2 = " -Kroupa -PARSEC -mcdata -mcseed="
	part3 = " -logterrsig="
	part4 = " -mbolerrsig="
	part5 = " > "+fullpath+"console"
	
	for i in range(0,50):
		rand = random.randint(0,5000)
		digit = '%02d' % (i,)
		outfile = open(fullpath+"mcseed_values.txt","a")#Storing the values of mcseed.
		out_string = "mcseed value " + str(digit) + "=" + str(rand) + "\n"
		outfile.write(out_string)
		outfile.close()
		comm = part1 + digit + part2 + str(rand) + part3 + str(lgsig) + part4 + str(mbol) + part5 + digit
		sp.call(comm.split())
	#combining the results using bestfit.
	#cmd10 = ["module","load","pgplot/5.2-sl6-gfortran"]
	#sp.call(cmd10)
	comarr = []
	comarr.append("zcombine "+fullpath+"out > "+fullpath+"out.zc")
	comarr.append("zcombine -unweighted -medbest -jeffreys -best="+fullpath+"out.zc "+fullpath+"out.mcmc > "+fullpath+"out.mcmc.zc")
	comarr.append("zcombine -unweighted -medbest -best="+fullpath+"out.zc "+fullpath+"out_?? > "+fullpath+"out.sys.zc")
	comarr.append("zcmerge "+fullpath+"out.zc "+fullpath+"out.mcmc.zc "+fullpath+"out.sys.zc -absolute > "+fullpath+"out.final")
	comarr.append("pg_cmd "+fullpath+"out.cmd "+fullpath+"out.ps -zclin="+fullpath+"out.final -zclog="+fullpath+"out.final")
	#combining wih hybrid uncertainities for comparison
	comarr.append("zcmerge "+fullpath+"out.zc "+fullpath+"out.mcmc.zc -absolute > "+fullpath+"out.hybrid.final")
	comarr.append("pg_cmd "+fullpath+"out.cmd "+fullpath+"out.hybrid.ps -zclin="+fullpath+"out.hybrid.final -zclog="+fullpath+"out.hybrid.final")
	for i in comarr:
		sp.call(i.split())
	
	#so we have created the pars file used in the main calcsfh runs, and completed a full calcsfh analysis of this galaxy
	#a plot has been created showing the fit and uncertainties

def fullFake(galdir, basis, pwd, galvals, goodfilt):
	#finds filter values that max. total lum. in output file. Uses this to find M/L ratio of galaxy
	
	#dumb paramater naming here
	#galdir -> basedir
	#basis -> scriptdir ???
	#pwd -> basedir/fakes
	#galvals -> flux, distance
	#goodfilter filter values to use
	
	#pars file split into three parts: header, filters, and footer
	#will first build these files to easily swap in new filter values as needed
	#first generate calcsfh pars file needed to create fake pars file
	sp.call(["mkdir",basis+"fakes"])
	makePars(galdir, pwd+"CparsBasis", goodfilt, "sfh_fullres")
	
	#use values to create start of fake pars file (up until timebins)
	p = open(pwd+"CparsBasis", 'r')	#open calcsfh pars file
	header = p.readline().split()	#read first line
	p.readline()			#read 2nd line, throw away
	head2 = p.readline().split()	#read 3rd line
	p.readline()			#don't care about lines 4,5
	p.readline()
	lunno = head2[0] #really don't know what this is
	dist = header[0] #grab dist and extint values from first line
	extin = header[3]
	p.close()
	g = open(pwd+"header", 'w')
	g.write(str(dist)+' '+str(extin)+' 0 '+str(lunno)+' -0.75'+"\n")
	g.write('2'+'\n')
	g.close()
	f = open(pwd+"filters", 'w')
	f.write(str(goodfilt[0])+' '+str(goodfilt[1])+' IRAC3.6'+'\n')
	f.write(str(goodfilt[2])+' '+str(goodfilt[3])+' IRAC4.5'+'\n')
	f.close()
	
	#open out.final from sfh_fullres, needed for sfr z values
	h = open(galdir+"sfh_fullres/out.hybrid.final","r")
	galmass = float(h.readline().split()[1])	#store gal mass for later, also sets file up to start reading timebins
	j = open("sfh_fullres", 'r') 	#open timebin file for # of timebins
	tbins = int(j.readline().split()[0])	#record number of timebins
	d = open(pwd+"footer", 'w')
	d.write(str(tbins)+'\n')	#write to fake pars file
	for x in range(0,tbins):	
		fLine = h.readline()	#read every output line in out.hybrid.final
		fStart = fLine[0:5] #grab start timebin
		fEnd = fLine[6:11] #grab end timebin
		fSFH = fLine[18:28] #grab SFH value
		fZee = fLine[51:57] #grab z value
		d.write("     "+fStart+" "+fEnd+" "+'%6.4f' % float(fSFH)+" "+'%4.1f' % float(fZee)+"\n")	#write to fake pars in required format
	
	#fake pars file now created. now we can close everything
	d.close()
	j.close()
	h.close()
		
	#run fake once to get starting lum value
	t = open(pwd+"parstest", 'w')
	with open(pwd+"header", 'r') as fobj:
		for line in fobj:
			t.write(line)
	with open(pwd+"filters", 'r') as fobj:
		for line in fobj:
			t.write(line)
	with open(pwd+"footer", 'r') as fobj:
		for line in fobj:
			t.write(line)
	comm = 'fake '+pwd+'parstest '+pwd+'outtest -full -KROUPA -PARSEC'
	sp.call(comm.split())
	sold = []
	sold.append(calclum(pwd+"outtest"))
	sold = sold + goodfilt	#list has starting lum and starting filter values
	
	#here is where we actually find the best filter values
	delta = .5		#choose how much values differ between runs
	runnum = 0		#keeps track of number of completed cycles
	while True:		#will stop loop 'manually' inside once end of cycle yield no positive change in lum
		if runnum > 100:
			print('you need a vacation')
			break
		flail = 0	#keeps track of permutation number
		lumlist = []	#records lum in given run
		permlist = []	#records filter values in given run	
		for w in range(-1,2):	#go through each perm of var inc/dec
			for x in range(-1,2):
				for y in range(-1,2):
					for z in range(-1,2):
						strflail = '%03d' % (flail,)		#convert run number to string for out files
						totest = sold[1:]			#grab values stored at end of previous cycle
						totest[0] = totest[0] + w * delta	#inc, dec, or stay constant depending on index values
						totest[1] = totest[1] + x * delta
						totest[2] = totest[2] + y * delta
						totest[3] = totest[3] + z * delta
						permlist.append(totest)
						makeFakePars(pwd, 'TEST'+strflail, totest)	#make pars file with 'test' to indicate temp file
						comm = 'fake '+pwd+'fakepars'+'TEST'+strflail+'.txt'+' '+pwd+'out'+'TEST'+strflail+' -full -KROUPA -PARSEC'
						sp.call(comm.split())
						flail = flail + 1
		for i in range(0, flail):
			strflail = '%03d' % (flail,)
			lumlist.append(calclum(pwd+'outTEST'+strflail, galvals[1]))	#record total luminosity in list entry
		maxloc, maxval = max(enumerate(lumlist), key=operator.itemgetter(1))		#find highest lum value and location after all trails complete
		runstr = '%03d' % (runnum,)			#convert cycle number to str for out files		
		if maxval > sold[0]:				#if highest found value, exceeds prev number, this is sucessful run
			sold[0] = maxval			#new stored lum value
			sold[1:] = permlist[maxloc]		#new stored filter values
			valstr = '%03d' % (maxloc,)		#store best run number as string for out files
			copyfile(pwd+'fakeparsTEST'+valstr+'.txt', pwd+'fakepars'+runstr+'.txt')	#copy temp pars of best run to new file
			copyfile(pwd+'outTEST'+valstr, pwd+'out'+runstr)				#copy temp out of best run to new file
			runnum = runnum + 1
		else:
			print('Best run found at index '+str(runnum - 1))			
			break
	
	#produce CMD of best run
	BestPlot(pwd, "out"+'%03d' % (runnum - 1,))
	
	#calculate mass/light ratio for galaxy
	#find total luminosity of best run
	totlum = calclum(pwd+"out"+'%03d' % (runnum - 1,), galvals[1])
	ratio = galmass/totlum	#In sol mass/sol lum
	
	#now compare magnitude from Spitzer with mag from fake luminosity
	SM, Sm = SpitMag(galvals)
	LM, Lm = LumMag(pwd+"out"+'%03d' % (runnum - 1,), galvals)

def SpitMag(galinfo):
	#finds magnitudes of gal from measured IRAC 3.6 flux
	fzero = 280.9	#Flux zero point in IRAC filter in Janskys
	fgal = galinfo[0]
	m = -2.5*math.log10(fgal/fzero)

	distMPC = galinfo[1]
	distPC = distMPC * 10 ** 6
	M = m - 5*(math.log10(distPC)-1)

	return M,m
	
def LumMag(path_to_fakeout, galinfo):
	#finds magnitudes of gal based on fake output
	lum = calclum(path_to_fakeout, galinfo[1])
	distMPC = galinfo[1]
	dist = distMPC * 10 ** 6
	logdist = 5*math.log10(dist/10)
	AbMag = 3.24 - 2.5*math.log10(lum)
	ApMag = AbMag + logdist
	return AbMag, ApMag

def BestPlot(pwd, outname):
	#produces CMD based on fake out file
	xVal = []	#IRAC3.6 - IRAC4.5
	yVal = []	#IRAC3.6
	with open(pwd+outname) as fobj:
		for line in fobj:
			row = line.split()		#read each line in fake out file
			yVal.append(float(row[0]))	#write each blue filter value to y-axis list
			comb = float(row[0])-float(row[1])	#find difference between red and blue mag on each line
			xVal.append(comb)		#write this value to x-axis lsit
	ax = plt.gca()		#allows for changing plot range
	ax.set_xlim([-0.2,0.05])	#focus x and y axis on main part of generated CMD. Will prob change with each galaxy
	ax.set_ylim([17,28])	
	ax.invert_yaxis()	#invert y-axis to match convention
	plt.scatter(xVal, yVal, s=0.2)
	plt.xlabel('IRAC3.6 - IRAC4.5')
	plt.ylabel('IRAC3.6')
	plt.title('CMD for best fake run: '+outname)
	plt.savefig(pwd+'CMD'+outname+'.png')	
	plt.close()
	
def makeFakePars(pwd, runber, totest):
	#creates pars file for fake program to use
	g = open(pwd+"fakepars"+runber, 'w')	#new pars file
	with open(pwd+"header", 'r') as fobj:
		for line in fobj:
			g.write(line)
	g.write(str(totest[0])+' '+str(totest[1])+' IRAC3.6'+'\n')	#write new filter values
	g.write(str(totest[2])+' '+str(totest[3])+' IRAC4.5'+'\n')
	with open(pwd+"footer", 'r') as fobj:
		for line in fobj:
			g.write(line)
	g.close()
	

def calclum(path_to_fakeout, galdist):
	#calculates lum of gal based on fake out file
	dist = galdist * 10 ** 6	#convert to Pcs
	logdist = 5*math.log10(dist/10)
	#Analyze test run to get starting luminosity
	mbol = []
	with open(path_to_fakeout) as fobj:
		for line in fobj:
			row = line.split()
			mbol.append(float(row[0]))	#only concerned with IRAC3.6 mag
	#convert mbol list to luminosity values
	for i in range(0, len(mbol) - 1):
		mbol[i] = 10 ** ((1/2.5)*(-mbol[i]+logdist+3.24))	#units of solar luminosity
	#sum luminosity
	tlum = 0
	for i in range(0, len(mbol) - 1):
		tlum = tlum + mbol[i]
	print(tlum)
	return tlum


def main():
	'''
	Edit the following:
	- GalName: Name of the folder of your galaxy in /acs/
	- GalPath: Path to the /galaxies/ dir
		I don't think you need it to be a path FROM your /oddraps/ folder to /galaxies/, but that is what I'm doing
	'''

	GalName = "10210_UGC8651"
	GalPath =  "../../../../../galaxies/"
	
	#find all cataloged infomation based on galaxy
	with open('GalCatalog','r') as fobj:
			for line in fobj:
				info = line.split()
				if info[0] == GalName:
					GalDir = info[1]
					GalFlux = float(info[2])
					GalDist = float(info[3])
					break
	basedir = GalPath+GalDir+"/"+GalName+"/metals_proc/"
	scriptr = basedir + "scriptdir/"
	fA, fB = setFolder(basedir)	#create work folder inside gal dir, move given files into it, send back filters used in file names
	print(fA,fB)
	#now we need to make any needed changes to these given files before using them in calcsfh
	Fstart = editFiles(scriptr, fA, fB)
	print(Fstart)
	#now run calcsfh with different filter values to find best depths
	#from here on out, we are running match commands and will need to use sbatch to run efficently
	bestDepth = calcFit(basedir, scriptr, Fstart)
	print(bestDepth)
	#now we can run the full calcsfh script for each timebin
	#fullCalc(scriptr, basedir+"sfh_fullres/", bestDepth, "sfh_fullres")
	#fullCalc(scriptr, basedir+"sfh_no_res/", bestDepth, "sfh_no_res")
	#fullCalc(scriptr, basedir+"sfh_starburst_v1res/", bestDepth, "sfh_starburst_v1res")
	#fullCalc(scriptr, basedir+"sfh_starburst_v2res/", bestDepth, "sfh_starburst_v2res")
	#all calcsfh runs have completed. Now to run fake to compute mass/light ratio
	#fullFake(basedir, scriptr, basedir+"fakes/", [GalFlux,GalDist], bestDepth)
	
if __name__ == "__main__":
    main()
