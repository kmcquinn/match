import time
import subprocess as sp
import sys
import random
import math
import operator
from shutil import copyfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fnmatch
import multiprocessing as mp
import os
import numpy as np
'''
HOW TO RUN ODDRAPS:
1.) Place lastest version in same folder as GalCatalog, sfh_full res file, teffdata file
2.) Edit sample batch script
	command: "python oddraps.py NAMEofGALAXYfolder -zinc=True/False"
	use your email, set time at ~48h
	may want to split calcFit and FullCalc into two separate jobs
		this is easily done by commenting out the functions you dont want for each script submission
	change job name, console output name accordingly
3.) Run by calling "sbatch SAMPLEscript"

Upon completion:
	/GALfolder/metals_proc/scriptdir/calctests has all cmds, ps files
	FiltResults file gives output name, filter depths, fit value
	/sfh_*timebin*/ folders have all fullCalc results.
		may need to run command to generate out.ps

oddraps = on dwarf disks, running a python script
'''
def setFolder(bir, phot, fake, pars):
	'''
	hey, this guy is taking the first file it sees that fits those wildcards
	that means if there is more than one match, any of those files is fair game for it to pull
	sorry
	'''
	#this creates a directory in the work folder and moves the given files (pars,phot,fake) over for safety
	#first setup folder in gal dir to build files
	sp.call(["mkdir",bir+"scriptdir/"])
	#find pars file in input_data folder
	if pars == None:
		fold = bir+"../conf_new_dol/"
		inList = sp.check_output(["ls",fold]).splitlines()
		cull = fnmatch.filter(inList,"*matchpars231")
		copyfile(fold+cull[0], bir+"scriptdir/basepars")
	else:
		copyfile(pars, bir+"scriptdir/basepars")
	
	fold = bir+"../proc_new_dol/"
	inList = sp.check_output(["ls",fold]).splitlines()
	
	if phot == None:
		exten = "*.gst.match"
		cull = fnmatch.filter(inList,exten)
		copyfile(fold+cull[0], bir+"scriptdir/phot")
	else:
		copyfile(phot, bir+"scriptdir/phot")
	if fake == None:
		cull = fnmatch.filter(inList,"*gst.matchfake")
		copyfile(fold+cull[0], bir+"scriptdir/fake")
	else:
		copyfile(fake, bir+"scriptdir/fake")
	#ok so all the files we need are in our new directory, and we grabbed the filter names from the filenames to compare with the given pars file. 

def editFiles(sir):
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

def makePars(basepath, newpath, depths, times, zinc):
	#produced calcsfh pars file based on given parameters
	#basepath: basedir
	#newpath: location, name of new pars file
	#depths: array of depth values used in this pars file
	#times: which timebin to append pars file with, assume this is a string with the name of the timebin file in /galaxies/ dir
	
	#open base pars file and write first 5 lines to new pars file
	g = open(basepath+'scriptdir/pars', 'r')
	p = open(newpath, 'w')
	p.write(g.readline())
	line2 = g.readline()[:-1]
	if zinc == True:
		line2 = "-2.0 0.1 0.1 -2.0 -0.9 -1.6 0.1"
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

def doWork(putin):
	flail, commdict = putin
	comm = commdict[flail][2]
	p = sp.Popen(comm.split(), stdout=sp.PIPE)
	out, err = p.communicate()
	while True:
		check = out.splitlines()[-1]
		if check.split()[0] == "Best":
			fit = float(check.split()[-1][4:])
			return fit
			break
		else:
			time.sleep(60)
			continue

def calcFit(bir, scr, filtStart, zinc):
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
	'''
	makePars(bir, scrstring+"parsFirst", sold[1:], "sfh_fullres")
	comm = "calcsfh "+scrstring+"parsFirst "+scrPhot+" "+scrFake+" "+scrstring+"outFirst -Kroupa -PARSEC"
	fitstr = sp.check_output(comm.split()).splitlines()[-1].split()[-1]
	#check consolefile to record starting fit value
	#format: 'fit=3573.515152'
	fitval = float(fitstr[4:])
	print("fitval is: "+str(fitval))
	#record this value to sold array
	sold[0] = fitval
	'''
	#let blue vary from 31 to 41
	#let red vary from 31 to 41
	
	delta = .33		#choose how much values differ between runs
	maxdelta = .25		#choose max deviation from start values
	runnum = 1		#keeps track of number of completed cycles
	maxrun = int(maxdelta/delta)
	
	flail = 0
	commdict = {} #commdict[flail] = [strflail, filtList, command, fit]
	runname = []
	fitlist = []
	g = open(scrstring+"FiltResults","w")
	g.write("Run Number\tDepth1\tDepth2\tFit Value\n")
	bluearr = np.arange(31.0,41.0,delta)
	redarr = np.arange(31.0,41.0,delta)
	for i in bluearr:
		for j in redarr:
			strflail = '%03d' % (flail,)		#convert run number to string for out files
			commdict[flail] = [strflail]
			parspath = scrstring+"calcparsTEST"+strflail	#paths to calcsfh files
			outpath = scrstring+"outTEST"+strflail
			consolepath = scrstring+"consoleTEST"+strflail
			totest = sold[1:]	#changing depth values up or down by delta
			totest[1] = i
			totest[3] = j
			commdict[flail].append(totest)
			makePars(bir, parspath, totest, "sfh_fullres", zinc)	#make pars file with 'test' to indicate temp file
			comm = "calcsfh "+parspath+" "+scrPhot+" "+scrFake+" "+outpath+" -Kroupa -PARSEC"	#command to send out
			if zinc == True:
				comm = comm + " -zinc"
			commdict[flail].append(comm)
			runname.append(flail)	#add runname to list
			flail = flail + 1
			#all commands created, now run them all with pool
	pool = mp.Pool(None)
	putin = [(i, commdict) for i in runname]
	result = pool.map_async(doWork, putin)
	pool.close()
	pool.join()
	outputfits = result.get()
	print(outputfits)
	print("length of outfits is "+str(len(outputfits)))
	print("length of runname is "+str(len(runname)))
	for i in range(0,len(outputfits)):
		commdict[i].append(outputfits[i])
		coolarr = commdict[i]
		outname = coolarr[0]
		Opath = scrstring+"outTEST"+outname
		sp.call(["pg_cmd",Opath+".cmd",Opath+".ps"])	#create output ps file for each run
		FiltA = str(coolarr[1][1])
		FiltB = str(coolarr[1][3])
		g.write(outname+"\t"+FiltA+"\t"+FiltB+"\t"+str(coolarr[3])+"\n")	#grab filter depths and write to results file
	minloc, minval = min(enumerate(outputfits), key=operator.itemgetter(1))
	beststr = '%03d' % (minloc,)
	g.write("Best run: outTEST"+beststr+", "+str(minval)+" with filter values "+str(commdict[minloc][1][1])+" "+str(commdict[minloc][1][3]))
	g.close()
	return commdict[minloc][1]
	
def Calcwork(arr):
	comm, output = arr[0], arr[1]
	f = open(output, "wb")	
	sp.call(comm.split(),stdout=f)
	f.close()
	return 0
def fullCalc(bpath, fullpath, goodDepths, tbins, zinc, mist):
	'''
	bpath -> /scriptdir/
	fullpath -> /sfh_fullres/
	'''
	comm = "mkdir "+fullpath
	sp.call(comm.split())
	#runs the full calcsfh workflow, incudes hybridMC, .ps plot of results
	photLoc = bpath+"phot"
	fakeLoc = bpath+"fake"
	#create needed pars file, run full calcsfh script
	makePars(bpath+"../", fullpath+"pars", goodDepths, tbins, zinc)
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
		fobj.readline()
		for line in fobj:
			row = line.split()
			vals.append([float(row[2]), float(row[3])])
			err.append(abs(Mv_diff - float(row[0])) + abs(Mi_diff - float(row[1])))
	minloc, minval = min(enumerate(err), key=operator.itemgetter(1))
	lgsig, mbol = vals[minloc][0], vals[minloc][1]
	#run calcsfh once for use with hybridMC
	comm1 = "calcsfh "+fullpath+"pars "+photLoc+" "+fakeLoc+" "+fullpath+"out -Kroupa "
	if mist == True:
		comm2 = "-MIST "
	else:
		comm2 = "-PARSEC "
	comm3 = "-mcdata"
	comm4 = ""
	if zinc == True:
		comm4 = " -zinc"
	comm = comm1 + comm2 + comm3 + comm4
	f = open(fullpath+"console.txt", "wb")	
	sp.call(comm.split(),stdout=f)
	f.close()
	#run hybridMC 1000 times
	comm = "hybridMC "+fullpath+"out.dat "+fullpath+"out.mcmc -tint=2.0 -nmc=10000 -dt=0.015"
	f = open(fullpath+"hybrid_console.txt", "wb")	
	sp.call(comm.split(),stdout=f)
	f.close()
	
	part1 = "calcsfh "+fullpath+"pars "+photLoc+" "+fakeLoc+" "+fullpath+"out_"
	part2 = " -Kroupa "
	if mist == True:
		part3 = "-MIST "
	else:
		part3 = "-PARSEC "
	if zinc == True:
		part4 = "-zinc -mcdata -mcseed="
	else:
		part4 = "-mcdata -mcseed="
	part5 = " -logterrsig="
	part6 = " -mbolerrsig="
	
	runarr = []
	for i in range(0,50):
		rand = random.randint(0,5000)
		digit = '%02d' % (i,)
		outfile = open(fullpath+"mcseed_values.txt","a")#Storing the values of mcseed.
		out_string = "mcseed value " + str(digit) + "=" + str(rand) + "\n"
		outfile.write(out_string)
		outfile.close()
		comm = part1 + digit + part2 + part3 + part4 + str(rand) + part5 + str(lgsig) + part6 + str(mbol)
		runarr.append([comm, fullpath+"console"+digit])
	pool = mp.Pool(None)
	pool.map_async(Calcwork, runarr)	#fills cpu cores with dowork jobs, each with different flail value from runname
	pool.close()
	pool.join()
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
		os.system(i)
	
	#so we have created the pars file used in the main calcsfh runs, and completed a full calcsfh analysis of this galaxy
	#a plot has been created showing the fit and uncertainties


def Fakework(comm):
	sp.call(comm.split())
	return 0
def fullFake(galdir, basis, pwd, galvals, goodfilt, zinc, mist, datdir):	
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
	galflux, galdist = galvals
	
	
	sp.call(["mkdir",pwd])
	makePars(galdir, pwd+"CparsBasis", goodfilt, "sfh_fullres", zinc)
	
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
	try:
		h = open(galdir+"sfh_fullresMIST/out.final","r")	
	except:
		h = open(galdir+"sfh_fullres/out.final","r")
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
	comm1 = 'fake '+pwd+'parstest '+datdir+'outtest -full -KROUPA '
	if mist == True:
		comm2 = "-MIST"
	else:
		comm2 = "-PARSEC"
	comm = comm1 + comm2
	sp.call(comm.split())
	sold = []
	sold.append(calclum(datdir+"outtest", galdist))
	sold = sold + goodfilt	#list has starting lum and starting filter values
	
	#here is where we actually find the best filter values
	delta = .1		#choose how much values differ between runs
	runnum = 0		#keeps track of number of completed cycles
	xdepth = []
	ydepth = []
	while True:		#will stop loop 'manually' inside once end of cycle yield no positive change in lum
		if runnum > 100:
			print('you need a vacation')
			print("Best run found at index 100")
			break
		flail = 0	#keeps track of permutation number
		lumlist = []	#records lum in given run
		permlist = []	#records filter values in given run
		commlist = []	
		for w in range(-50,50):	#go through each perm of var inc/dec
			#for x in range(-12,15):
			strflail = '%03d' % (flail,)		#convert run number to string for out files
			totest = sold[1:]			#grab values stored at end of previous cycle
			totest[1] = totest[1] + w * delta	#inc, dec, or stay constant depending on index values
			#totest[3] = totest[3] + x * delta
			totest[3] = 34.0
			permlist.append(totest)
			makeFakePars(pwd, 'TEST'+strflail, totest)	#make pars file with 'test' to indicate temp file
			comm1 = 'fake '+pwd+'fakepars'+'TEST'+strflail+' '+datdir+'out'+'TEST'+strflail+' -full -KROUPA '
			if mist == True:
				comm2 = "-MIST"
			else:
				comm2 = "-PARSEC"
			comm = comm1 + comm2
			commlist.append(comm)
			flail = flail + 1
		pool = mp.Pool(None)
		pool.map_async(Fakework, commlist)
		pool.close()
		pool.join()	
		for i in range(0, flail):
			strflail = '%03d' % (i,)
			lum = calclum(datdir+'outTEST'+strflail, galdist)
			rfilt = permlist[i][1]
			xdepth.append(rfilt)
			ydepth.append(lum)
			lumlist.append(lum)	#record total luminosity in list entry
		maxloc, maxval = max(enumerate(lumlist), key=operator.itemgetter(1))		#find highest lum value and location after all trails complete
		runstr = '%03d' % (runnum,)			#convert cycle number to str for out files		
		if maxval > sold[0]:				#if highest found value, exceeds prev number, this is sucessful run
			sold[0] = maxval			#new stored lum value
			sold[1:] = permlist[maxloc]		#new stored filter values
			valstr = '%03d' % (maxloc,)		#store best run number as string for out files
			copyfile(pwd+'fakeparsTEST'+valstr, pwd+'fakepars'+runstr)	#copy temp pars of best run to new file
			copyfile(datdir+'outTEST'+valstr, pwd+'out'+runstr)				#copy temp out of best run to new file
			runnum = runnum + 1
		else:
			print('Best run found at index '+str(runnum - 1))			
			break
	
	#produce CMD of best run
	g = open(pwd+"results","a")
	PlotCurve(pwd,xdepth,ydepth)
	BestPlot(pwd, "out"+'%03d' % (runnum - 1,))
	g.write("CMD of best run created at out"+'%03d' % (runnum - 1,)+"\n")
	g.write("curve of growth plotted at CMDcurve.png"+"\n")
	#calculate mass/light ratio for galaxy
	#find total luminosity of best run
	totlum = calclum(pwd+"out"+'%03d' % (runnum - 1,), galdist)
	g.write("Total Luminosity is "+str(totlum)+"\n")
	g.write("Gal Mass is "+str(galmass)+"\n")
	ratio = galmass/totlum	#In sol mass/sol lum
	g.write("Mass/Light ratio is "+str(ratio)+"\n")
	#now compare magnitude from Spitzer with mag from fake luminosity
	SM, Sm = SpitMag(galvals)
	g.write("Spitzer absol and appar mag: "+str(SM)+" "+str(Sm)+"\n")
	LM, Lm = LumMag(pwd+"out"+'%03d' % (runnum - 1,), galvals)
	g.write("Fake absol and appar mag: "+str(LM)+" "+str(Lm)+"\n")
	
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
	
def PlotCurve(pwd, xdepth, ydepth):
	plt.scatter(xdepth, ydepth, s=0.2)
	plt.xlabel('Blue Filter Depth')
	plt.ylabel('Total Luminosity')
	plt.title('Curve of growth w/ varying depth')
	plt.savefig(pwd+'CMDcurve.png')	
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
			try:
				if float(row[0]) < 98.0 and float(row[1]) < 98.0:
					mbol.append(float(row[0]))	#only concerned with IRAC3.6 mag
			except:
				print("for some reason this happened:", row)
	#convert mbol list to luminosity values
	for i in range(0, len(mbol)):
		mbol[i] = 10 ** ((1/2.5)*(-mbol[i]+logdist+3.24))	#units of solar luminosity
	#sum luminosity
	tlum = 0
	for i in range(0, len(mbol)):
		tlum = tlum + mbol[i]
	print(tlum)
	return tlum


def main():
	'''
	Your galaxy run needs to change:
	1.) GalName to the folder name of your galaxy
	2.) basedir: change if you happen to be in the /oddraps/ folder in a different directory that what I assume.
		-Change so current directory connects to metals_proc for your gal as shown
	3.) Comment out or add what functions needed for your analysis
		-setFolder: creates folder for our use, moves given phot, fake, pars files. Returns filter names given in phot file name
		-editFiles: edits given pars file to work with current calcsfh version. overwrites filter names in pars file if different to file name. Returns filter depths found in given pars file
		-calcFit: Generates a bunch of different calcsfh runs with different filter depths centered around those given. Produces FiltResults file to connect output files with filter and fit values. Also finds the run with the lowest fit value. Returns the filter values that produced the lowest fit
		-fullCalc: MC runs for given filter depths and timebin. Also combines and produces all needed .ps outputs
		-fullFake: Finds IRAC filter depths that produce highest "fake" luminosity for galaxy. Compares with known magnitude of gal, produces M/L ratio
	    You will need setFolder and editFiles the first time on a given galaxy in order to run any other functions
	    Information from each function is not ingelligently stored (ie saved in reference file to save time on repeat runs). I'll add this once everything works.
	'''
	#python oddraps.py GalFolder -zinc=True/False -time=full/no/v1/v2 -phot=PhotPath -fake=FakePath -pars=ParsPath -mist=True/False
	#pars for fullCalc: fullCalc(scriptr, basedir+tibin+"/", bestDepth, tibin=", Zinc, Mist)
	'''
	adding mist to oddraps:
		2. Change fullCalc/fullFake to change comm's based on zinc flag
	'''
	Zinc = False
	tibin = "sfh_fullres"
	phot = None
	fake = None
	pars = None
	mist = False
	for i in range(1,len(sys.argv)):
		if i == 1:
			GalName = sys.argv[i]
		if sys.argv[i][1:5] == "zinc":
			Zinc = sys.argv[i][6:]
			if Zinc == "True":
				Zinc = True
			else:
				Zinc = False
		if sys.argv[i][1:5] == "time":
			tibin = sys.argv[i][6:]	#form: -time=full/no/v1/v2
			if tibin == "full":
				tibin = "sfh_fullres"
			elif tibin == "no":
				tibin = "sfh_no_res"
			elif tibin == "v1":
				tibin = "sfh_starburst_v1res"
			elif tibin == "v2":
				tibin = "sfh_starburst_v2res"
		if sys.argv[i][1:5] == "phot":
			phot = sys.argv[i][6:]
		if sys.argv[i][1:5] == "fake":
			fake = sys.argv[i][6:]
		if sys.argv[i][1:5] == "pars":
			pars = sys.argv[i][6:]
		if sys.argv[i][1:5] == "mist":
			mist = sys.argv[i][6:]
			if mist == "True":
				mist = True
			else:
				mist = False
	#find all cataloged infomation based on galaxy
	runFit = 1
	with open('GalCatalog','r') as fobj:
			for line in fobj:
				info = line.split()
				if info[0] == GalName:
					GalDir = info[1]
					GalFlux = float(info[2])
					GalDist = float(info[3])
					if len(info) > 4:
						bestDepth = info[4:8]
						bestDepth = [float(i) for i in bestDepth]
						runFit = 0
					break
	basedir = "/work/04316/kmcquinn/wrangler/metals/galaxies/"+GalDir+"/"+GalName+"/metals_proc/"
	scriptr = basedir + "scriptdir/"
	setFolder(basedir, phot, fake, pars)	#create work folder inside gal dir, move given files into it, send back filters used in file names
	#now we need to make any needed changes to these given files before using them in calcsfh
	Fstart = editFiles(scriptr)
	#print(Fstart)
	#now run calcsfh with different filter values to find best depths
	#from here on out, we are running match commands and will need to use sbatch to run efficently
	if runFit == 1:
		bestDepth = calcFit(basedir, scriptr, Fstart, Zinc)
	#print(bestDepth)
	#now we can run the full calcsfh script for each timebin
	if mist == True:
		calcfolder = basedir+tibin+"MIST"+"/"
		fakefolder = basedir+"fakesMIST/"
	else:
		calcfolder = basedir+tibin+"/"
		fakefolder = basedir+"fakes/"
	#fullCalc(scriptr, calcfolder, bestDepth, tibin, Zinc, mist)
	datName = GalName+"_store"
	comm = "mkdir $DATA/"+datName
	datdir = "$DATA/"+GalName+"_store/"
	sp.call(comm.split())
	fullFake(basedir, scriptr, fakefolder, [GalFlux,GalDist], bestDepth, Zinc, mist, datdir)
	
if __name__ == "__main__":
    main()
