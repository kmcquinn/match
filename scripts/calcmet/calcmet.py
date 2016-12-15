
# coding: utf-8

# In[38]:

import numpy as np
from sys import argv as arg
import os
from astropy.table import Table,Column,hstack,vstack
zsol = 8.69     # Solar z
fake = 8.0     # Test gas z


# In[39]:

# Searches for item in a list and returns the index of the item in the list if it is in the list
def itin(x,pop):
    for i in range(len(pop)):
        if x == pop[i]:
            return int(i)

'''
Defines variables from out.final file in the table 'dat'.
'''

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
    
'''
Reads in oxygen abundances from galaxy_list_O.txt as the table 'gas'
'''

def gastab(boop):
    global gas
    f=open(boop,'r')
    nam=f.readline().split()
    gal=[[i.split()[0],i.split()[3],i.split()[5],i.split()[6],i.split()[8]] for i in f]
    f.close()
    gas=Table(rows=gal,names=['gal','z','zerr','NO','NOerr'],dtype=['S','f8','f8','f8','f8'])
    return gas

'''
Reads in HI abundances from metals_opticaldata.txt as the table 'hi'
'''

def hitab(boop):
    global hi
    f=open(boop,'r')
    for i in range(31):
        f.readline()
    ale=[[i.split()[1], i.split()[18], i.split()[19]] for i in f]     # Takes galaxy name, 
                                                                      # log HI flux, and dist in Mpc
    f.close()
    hi=Table(rows=ale,names=['gal','loghi','dist'],dtype=['S','f8','f8'])


# In[40]:

'''
Calculates the amount of oxygen formed in the galaxy with 3 different nucleosynthesis yields from 
Nomoto et al. 2006, which are p1, p2, and p3.
'''

def cogal(t,u,l):     # t = total star formation
    t=float(t)
    p1=0.0054
    p2=0.00658
    p3=0.0086
    return [[t*p1,u*p1,l*p1],[t*p2,u*p2,l*p2],[t*p3,u*p3,l*p3]]


# In[41]:

'''
Calculates the mass of oxygen in the gas.
'''
# Need to include zerr in final output!

def cogas(nam):     # nam = galaxy name, no 0's before number (ex. UGC8508, not UGC08508)

    if len(nam)!= 3 and nam[3] == ' ':
        nam = nam[0:3]+nam[4:]
############### Takes O abundance from known data ###################

    if nam in gas['gal']:     # Checks if galaxy's oxygen abundances are available
        i=itin(nam,gas['gal'])
        z=gas['z'][i]
        zerr=gas['zerr'][i]

    else:     # Changes input nam if galaxy isn't found in data
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

############### Calculates atomic hydrogen gas mass from HI flux ##################

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
    return [gm*16.*10**(z-12.),abs(gm*16.*10**(z-12.)*np.log(10))*zerr]


# In[42]:

'''
Calculates the amount of oxygen locked in the stars.
'''

def costar(r,ru,rl,st,et,z,zerru,zerrl):    # Takes sfr, sfr upper and lower uncertainty, time bin, metallicity,
                                            # upper and lower uncertainties
    Rec = 1./3.     # Recycling Fraction
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


# In[43]:

'''
Calculates the total oxygen budget.
'''
def obud(g,s,t):     # Takes oxygen in gas, stars, and total oxygen formed
    ans=[]
    erru=[]
    errl=[]
    for i in range(len(t)):
        ans.append((g[0]+s[0])/t[i][0])
        erru.append((((g[1]**2+s[1]**2)**(1./2.)/(g[0]+s[0]))**2+(t[i][1]/t[i][0])**2)**(1./2.)*ans[i])
        errl.append((((g[1]**2+s[1]**2)**(1./2.)/(g[0]+s[0]))**2+(t[i][1]/t[i][0])**2)**(1./2.)*ans[i])
    return ans,erru,errl


# In[44]:

'''
Add results to a really ugly file.
'''
# apparently blank.ljust is a good way of making a really pretty file
def maketab(nam,filnam,output):
    define(filnam)
    gastab('galaxy_list_O.txt')
    hitab('metals_opticaldata.txt')
    a=cogas(nam)
    b=costar(dat['sfr'],dat['sfru'],dat['sfrl'],dat['start'],dat['end'],dat['met'],dat['metu'],dat['metl'])
    c=cogal(totsf,totsfu,totsfl)
    d=obud(a,b,c)
    if a!='No oxygen abundance available for '+nam+'.' and a!='No HI flux available for '+nam+'.':
        if output in os.listdir("."):
            thi=open(output,'a')
            thi.write('\n')
            thi.write(nam+'\t'+str(d[0][1])+'\t'+str(d[1][1])+'\t'+str(d[2][1])+'\t'+str(a[0])+'\t'+str(a[1])+'\t'+str(b[0])+'\t'+str(b[1])
                      +'\t'+str(b[2])+str(c[0][0])+'\t'+str(c[0][1])+'\t'+str(c[0][2])+'\t'+str(c[1][0])+'\t'
                      +str(c[1][1])+'\t'+str(c[1][2])+'\t'+str(c[2][0])+'\t'+str(c[2][1])+'\t'+str(c[2][2]))
            thi.close()
        else:
            thi=open(output,'w')
            thi.write('Name\tMRF\t+\t-\tStellar Mass\t+\t-\tO_gas\t+/-\tO_star\t+\t-\tTot1\t+\t-\tTot2\t+\t-\tTot3\t+\t-')
            thi.write('\n')
            thi.write(nam+'\t'+str(d[0][1])+'\t'+str(d[1][1])+'\t'+str(d[2][1])+str(totsf)+'\t'+str(totsfu)+'\t'+str(totsfl)+'\t'+str(a[0])+'\t'+str(a[1])+'\t'+str(b[0])+'\t'+str(b[1])
                      +'\t'+str(b[2])+'\t'+str(c[0][0])+'\t'+str(c[0][1])+'\t'+str(c[0][2])+'\t'+str(c[1][0])+'\t'
                      +str(c[1][1])+'\t'+str(c[1][2])+'\t'+str(c[2][0])+'\t'+str(c[2][1])+'\t'+str(c[2][2]))
            thi.close()


# In[45]:

# Proper syntax: python calcmet.py [name of file with data in it]

def main():
        res="sfh_fullres"
        galdir='/work/04316/kmcquinn/wrangler/metals/galaxies/acs'
        nam=open('filnam','r')
        with open('filnam') as f:
                bu=f.read().splitlines()
        naml=[]
        for i in bu:
                naml.append(i.split("\t"))
        
        mist= "False"
        for g in arg:
            if g[:6]=='-mist=':
                mist=g[6:]
        if mist == "True":
            for i in naml:
                maketab(i[1],galdir+'/'+i[0]+'/metals_proc/'+res+'MIST/out.final','mist'+res+'_all')
                maketab(i[1],galdir+'/'+i[0]+'/metals_proc/'+res+'MIST/out.hybrid.final','mist'+res+'_nosys')
        if mist == "False":
            for i in naml:
                maketab(i[1],galdir+'/'+i[0]+'/metals_proc/'+res+'/out.final','parsec'+res+'_all')
                maketab(i[1],galdir+'/'+i[0]+'/metals_proc/'+res+'/out.hybrid.final','parsec'+res+'_nosys')
    

if __name__ == '__main__':
    main()


# In[ ]:



