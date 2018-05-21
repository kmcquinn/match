import numpy as np
from sys import argv as arg
import os
import argparse
from astropy.table import Table,Column,hstack,vstack
zsol = 8.69     # Solar z
fake = 8.0     # Test gas z

'''
SYNTAX:
python calcmet.py [file with galaxies] [resolution] [library]
'''

'''GENERAL FUNCTIONS'''
# Searches for item in a list and returns the index of the item in the list if it is in the list
def itin(x,pop):
    for i in range(len(pop)):
        if x == pop[i]:
            return int(i)

def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('galaxy_list',action='store',help='file containing list of galaxies whose MRFs will be calculated')
    parser.add_argument('resolution',action='store',help='resolution to be used for MRF calculation')
    parser.add_argument('library',action='store',help='stellar evolution library used in desired run')
    parser.add_argument('--camera','-c',dest='camera',action='store',help='camera used for observations',default='acs')
    args = parser.parse_args()
    return args


'''LOADING DATA'''
#Defines variables from out.final file in the table 'dat'.
def define(bop):
    global totsf
    global totsfu
    global totsfl
    global dat
    
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
    
#Reads in oxygen abundances from galaxy_list_O.txt as the table 'gas'
def gas_table(boop):
    global gas
    f=open(boop,'r')
    nam=f.readline().split()
    gal=[[i.split()[0],i.split()[3],i.split()[5],i.split()[6],i.split()[8]] for i in f]
    f.close()
    gas=Table(rows=gal,names=['gal','z','zerr','NO','NOerr'],dtype=['S','f8','f8','f8','f8'])
    return gas

#Reads in HI abundances from metals_opticaldata.txt as the table 'hi'
def HI_table(boop):
    global hi
    f=open(boop,'r')
    for i in range(31):
        f.readline()
    ale=[[i.split()[1], i.split()[18], i.split()[19]] for i in f]     # Takes galaxy name, 
                                                                      # log HI flux, and dist in Mpc
    f.close()
    hi=Table(rows=ale,names=['gal','loghi','dist'],dtype=['S','f8','f8'])


'''MRF CALCULATIONS'''
# Calculates the amount of oxygen formed in the galaxy with 3 different nucleosynthesis yields from 
#    Nomoto et al. 2006, which are p1, p2, and p3.
def Omass_total(t,u,l):     # t = total star formation
    t=float(t)
    p1=0.0054
    p2=0.00658
    p3=0.0086
    return [[t*p1,u*p1,l*p1],[t*p2,u*p2,l*p2],[t*p3,u*p3,l*p3]]

# Calculates the mass of oxygen in the gas.
# Need to include zerr in final output!
def Omass_gas(nam):     # nam = galaxy name
    if len(nam)!= 3 and nam[3] == ' ':
        nam = nam[0:3]+nam[4:]

    # Takes O abundance from known data
    if nam in gas['gal']:     # Checks if galaxy's oxygen abundances are available
        i=itin(nam,gas['gal'])
        z=gas['z'][i]
        zerr=gas['zerr'][i]
    else:     # Changes input name if galaxy isn't found in data
        mam = nam[0:3]+'0'+nam[3:]
        i=itin(mam,gas['gal'])
        if i != None:
            z=gas['z'][i]
            zerr=gas['zerr'][i]
        else:
            mam = mam[0:3]+'0'+mam[3:]
            i=itin(mam,gas['gal'])
            if i != None:
                z=gas['z'][i]
                zerr=gas['zerr'][i]
            else:
                return 'No oxygen abundance available for '+nam+'.'

    # Calculates atomic hydrogen gas mass from HI flux
    if nam in hi['gal']:
        i=itin(nam,hi['gal'])
        ahg=2.356e5*(float(hi['dist'][i])**2)*10**(float(hi['loghi'][i]))
        # ahg = atomic hydrogen gas mass
    else:     # Changes input nam if galaxy isn't found in database
        ham = nam[0:3]+'0'+nam[3:]
        i=itin(ham,hi['gal'])
        if i!=None:
            ahg=2.356e5*(float(hi['dist'][i])**2)*10**(float(hi['loghi'][i]))
        else:
            ham = ham[0:3]+'0'+ham[3:]
            i=itin(ham,gas['gal'])
	    if i != None:
                ahg=2.356e5*(float(hi['dist'][i])**2)*10**(float(hi['loghi'][i]))
            else:
                return 'No HI flux available for '+nam+'.'

    agm = 1.33 * ahg     # Total atomic gas mass, includes helium
    mg = 0.1 * agm     # Molecular gas mass (assuming no availabe measurements)
    gm = agm + mg     # Total mass of gas
    return [[gm*16.*10**(z-12.),abs(gm*16.*10**(z-12.)*np.log(10))*zerr],ahg]

#Calculates the amount of oxygen locked in the stars.
def Omass_stars(r,ru,rl,st,et,z,zerru,zerrl):    # Takes sfr, sfr upper and lower uncertainty, time bin, metallicity, upper and lower uncertainties
    Rec = .436    # Recycling Fraction
    somd = []     # Scaled oxygen mass density for each time bin
    somderu = []
    somderl = []
    sm = []     # Stellar mass for each time bin
    smeru = []
    smerl = []
    for i in range(len(z)):
        if r[i] != 0 and z[i] != 0:
            ond = z[i]+zsol     # Oxygen number density in stars
            osca = 1.     # Oxygen scale factor (need to find actual value)
            somd.append(10**(ond-12+np.log10((16.)/(0.75*1.0079 + 0.25*4.0026)))*osca)
            # Next 2 lines are error propagation
            somderu.append((10**(ond-12+np.log10((16.)/(0.75*1.0079 + 0.25*4.0026)))*np.log(10)*osca)*zerru[i])
            somderl.append((10**(ond-12+np.log10((16.)/(0.75*1.0079 + 0.25*4.0026)))*np.log(10)*osca)*zerrl[i])
            sm.append((10**et[i]-10**st[i])*r[i])
            # Next 2 lines are error propagation
            smeru.append((10**et[i]-10**st[i])*ru[i])
            smerl.append((10**et[i]-10**st[i])*rl[i])
    om = []     # Oxygen mass in stars for each time bin
    omerru = []
    omerrl = []
    for i in range(len(somd)):
        #Next 2 lines are error propagation
        om.append((1-Rec)*sm[i]*somd[i])
        omerru.append(om[i]*(1-Rec)*((smeru[i]/sm[i])**2.+(somderu[i]/somd[i])**2.)**(1./2.))
        omerrl.append(om[i]*(1-Rec)*((smerl[i]/sm[i])**2.+(somderl[i]/somd[i])**2.)**(1./2.))
    return [sum(om), (sum([i**2 for i in omerru]))**(1./2.), (sum([i**2 for i in omerrl]))**(1./2.)]     # Total oxygen mass 
                                                                                                         # in stars

# Calculates the total oxygen budget.
def calculate_mrf(g,s,t):     # Takes oxygen in gas, stars, and total oxygen formed
    tans=[] ; gans=[] ; sans=[]
    terru=[] ; gerru=[] ; serru=[]
    terrl=[] ; gerrl=[] ; serrl=[]
    # add in factor of 2 uncertainty in stellar o abundance upper uncert based on O/Fe ratios McWilliam 1997 AARA
    s[1] = (s[1]**2 + 0.5**2)**(1./2.)
    for i in range(len(t)):
        tans.append((g[0]+s[0])/t[i][0])
        gans.append(g[0]/t[i][0])
        sans.append(s[0]/t[i][0])
        terru.append((((g[1]**2+s[1]**2)**(1./2.)/(g[0]+s[0]))**2+(t[i][1]/t[i][0])**2)**(1./2.)*tans[i])
        terrl.append((((g[1]**2+s[2]**2)**(1./2.)/(g[0]+s[0]))**2+(t[i][2]/t[i][0])**2)**(1./2.)*tans[i])
        gerru.append(((g[1]/(g[0]))**2+(t[i][1]/t[i][0])**2)**(1./2.)*gans[i])
        gerrl.append(((g[1]/(g[0]))**2+(t[i][2]/t[i][0])**2)**(1./2.)*gans[i])
        serru.append(((s[1]/(s[0]))**2+(t[i][1]/t[i][0])**2)**(1./2.)*sans[i])
        serrl.append(((s[2]/(s[0]))**2+(t[i][2]/t[i][0])**2)**(1./2.)*sans[i])
    return tans,terru,terrl,gans,gerru,gerrl,sans,serru,serrl


'''SAVE RESULTS'''
# Add results to a file.
# apparently blank.ljust is a good way of making a really pretty file
def save_data(nam,filnam,output):
    define(filnam)
    gas_table('galaxy_list_O.txt')
    HI_table('metals_opticaldata.txt')
    a=Omass_gas(nam)[0]
    ahg = Omass_gas(nam)[1]
    b=Omass_stars(dat['sfr'],dat['sfru'],dat['sfrl'],dat['start'],dat['end'],dat['met'],dat['metu'],dat['metl'])
    c=Omass_total(totsf,totsfu,totsfl)
    d=calculate_mrf(a,b,c)
    if a!='No oxygen abundance available for '+nam+'.' and a!='No HI flux available for '+nam+'.':
        if output in os.listdir("."):
            thi=open(output,'a')
            thi.write('\n')
            thi.write(nam+'\t'+str(d[0][1]*100)+'\t'+str(d[1][1]*100)+'\t'+str(d[2][1]*100)+'\t'+str(d[3][1]*100)+'\t'+str(d[4][1]*100)+'\t'+str(d[5][1]*100)+'\t'+str(d[6][1]*100)+'\t'+str(d[7][1]*100)+'\t'+str(d[8][1]*100)+'\t'+str(totsf)+'\t'+str(ahg)+'\t'+str(ahg/totsf))
            thi.close()
        else:
            thi=open(output,'w')
            thi.write('GalaxyName\tMRF\t+\t-\tOgas\t+\t-\tOstars\t+\t-\tMstar\tMgas\tMgas/Mstar')
            thi.write('\n')
            thi.write(nam+'\t'+str(d[0][1]*100)+'\t'+str(d[1][1]*100)+'\t'+str(d[2][1]*100)+'\t'+str(d[3][1]*100)+'\t'+str(d[4][1]*100)+'\t'+str(d[5][1]*100)+'\t'+str(d[6][1]*100)+'\t'+str(d[7][1]*100)+'\t'+str(d[8][1]*100)+'\t'+str(totsf)+'\t'+str(ahg)+'\t'+str(ahg/totsf))
	    thi.close()


# Proper syntax: python calcmet.py [name of file with data in it]
def main():
    # Parses through command line arguments
    args = parse_options()
    galdir='/work/04316/kmcquinn/wrangler/metals/galaxies/'+args.camera
    nam=open(args.galaxy_list,'r')
    with open(args.galaxy_list) as f:
            bu=f.read().splitlines()
    naml=[]
    for i in bu:
            naml.append(i.split("\t"))
    for i in naml:
            save_data(i[1],galdir+'/'+i[0]+'/metals_proc/'+args.resolution+'_'+args.library+'/out.hybrid.final',args.library+'_mrflist')



if __name__ == '__main__':
    main()
