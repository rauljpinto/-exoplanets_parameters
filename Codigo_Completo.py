#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import astropy
from astropy.table import Table
from astropy import constants as const

t = Table.read(r'C:\Users\Raul\Desktop\Exoplanets\Complete_table.csv', format='ascii.csv')

sb = const.sigma_sb
pi = np.pi
Ls = const.L_sun
au = 149597870700

R_est = t['st_rad']*695500000
a = t['pl_orbsmax']*au
dist = t['st_dist']*(3.086*(10**16))
AngDiam = (2*R_est)/dist
Met = t['st_metfe']


# In[2]:


# (B - V) Band

B_V = (t['st_bj']-t['st_vj'])
BV_eff = 0.5725 + 0.4722*B_V + 0.0086*(B_V**2) - 0.0628*B_V*Met - 0.0038*Met - 0.0051*(Met**2)
Teff1 = 5040/BV_eff

#Equilibrium Temperature
Teq = (Teff1*0.9147)*(R_est/(2.0*a))**(1.0/2.0)

#Luminosity
L = 4*pi*(R_est**2.0)*sb*(Teff1)**4.0

#Insolation Flux
S = ((L/Ls)*(au/a)**2.0)*1.361

#Habitable Zone
HZ_in = (0.75*((L/Ls)**0.5))*au
HZ_center = (1.0*((L/Ls)**0.5))*au
HZ_out = (1.77*((L/Ls)**0.5))*au
deltaHZ = HZ_out - HZ_in


fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Equilibrium Temperature (K)')
ax.grid()

ax.scatter(Teff1, Teq, s=1)
plt.xlim(1500, 8000)
plt.ylim(0, 2500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Insolation Flux (Earth Flux)')
ax.grid()

ax.scatter(Teff1, S, s =1)
plt.xlim(2000, 7000)
plt.ylim(-10, 500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Habitable Zone (m)')
ax.grid()

ax.scatter(Teff1, HZ_center,  s = 1)
plt.xlim(0, 8000)
plt.ylim(0, 1e12)
plt.show()
plt.clf()


# In[3]:


# (V - J) Band

V_J = (t['st_vj']-t['st_j'])
VJ_eff = 0.4997 + 0.3504*V_J - 0.0230*(V_J**2) - 0.0295*V_J*Met + 0.0468*Met + 0.0037*(Met**2)
Teff4 = 5040/VJ_eff

#Equilibrium Temperature
Teq = (Teff4*0.9147)*(R_est/(2.0*a))**(1.0/2.0)

#Luminosity
L = 4*pi*(R_est**2.0)*sb*(Teff4)**4.0

#Insolation Flux
S = ((L/Ls)*(au/a)**2.0)*1.361

#Habitable Zone
HZ_in = (0.75*((L/Ls)**0.5))*au
HZ_center = (1.0*((L/Ls)**0.5))*au
HZ_out = (1.77*((L/Ls)**0.5))*au
deltaHZ = HZ_out - HZ_in


fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Equilibrium Temperature (K)')
ax.grid()

ax.scatter(Teff4, Teq, s=1)
plt.xlim(1500, 8000)
plt.ylim(0, 2500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Insolation Flux (Earth Flux)')
ax.grid()

ax.scatter(Teff4, S, s =1)
plt.xlim(2000, 7000)
plt.ylim(-10, 500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Habitable Zone (m)')
ax.grid()

ax.scatter(Teff4, HZ_center,  s = 1)
plt.xlim(0, 8000)
plt.ylim(0, 1e12)
plt.show()
plt.clf()


# In[4]:


# (V - H) Band

V_H = (t['st_vj']-t['st_h'])
VH_eff = 0.5341 + 0.2517*V_H - 0.0100*(V_H**2) - 0.0236*V_H*Met + 0.0523*Met + 0.0044*(Met**2)
Teff5 = 5040/VH_eff

#Equilibrium Temperature
Teq = (Teff5*0.9147)*(R_est/(2.0*a))**(1.0/2.0)

#Luminosity
L = 4*pi*(R_est**2.0)*sb*(Teff5)**4.0

#Insolation Flux
S = ((L/Ls)*(au/a)**2.0)*1.361

#Habitable Zone
HZ_in = (0.75*((L/Ls)**0.5))*au
HZ_center = (1.0*((L/Ls)**0.5))*au
HZ_out = (1.77*((L/Ls)**0.5))*au
deltaHZ = HZ_out - HZ_in


fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Equilibrium Temperature (K)')
ax.grid()

ax.scatter(Teff5, Teq, s=1)
plt.xlim(1500, 8000)
plt.ylim(0, 2500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Insolation Flux (Earth Flux)')
ax.grid()

ax.scatter(Teff5, S, s =1)
plt.xlim(2000, 7000)
plt.ylim(-10, 500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Habitable Zone (m)')
ax.grid()

ax.scatter(Teff5, HZ_center,  s = 1)
plt.xlim(0, 8000)
plt.ylim(0, 1e12)
plt.show()
plt.clf()


# In[6]:


# (J - Ks) Band

J_Ks = (t['st_j']-t['st_k'])
JKs_eff = 0.6524 + 0.5813*J_Ks + 0.1225*(J_Ks**2) - 0.0646*J_Ks*Met + 0.0370*Met + 0.0016*(Met**2)
Teff7 = 5040/JKs_eff

#Equilibrium Temperature
Teq = (Teff7*0.9147)*(R_est/(2.0*a))**(1.0/2.0)

#Luminosity
L = 4*pi*(R_est**2.0)*sb*(Teff7)**4.0

#Insolation Flux
S = ((L/Ls)*(au/a)**2.0)*1.361

#Habitable Zone
HZ_in = (0.75*((L/Ls)**0.5))*au
HZ_center = (1.0*((L/Ls)**0.5))*au
HZ_out = (1.77*((L/Ls)**0.5))*au
deltaHZ = HZ_out - HZ_in


fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Equilibrium Temperature (K)')
ax.grid()

ax.scatter(Teff7, Teq, s=1)
plt.xlim(1500, 8000)
plt.ylim(0, 2500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Insolation Flux (Earth Flux)')
ax.grid()

ax.scatter(Teff7, S, s =1)
plt.xlim(2000, 7000)
plt.ylim(-10, 500)
plt.show()
plt.clf()

fig, ax = plt.subplots()

ax.set(xlabel='Star Effective Temperature (K)', ylabel='Habitable Zone (m)')
ax.grid()

ax.scatter(Teff7, HZ_center,  s = 1)
plt.xlim(0, 8000)
plt.ylim(0, 1e12)
plt.show()
plt.clf()


# In[6]:


JK = np.array(Teff7)
JKlst = JK.tolist()
VJ = np.array(Teff4)
VJlst = VJ.tolist()


ListaJK=[]
ListaVJ=[]


for i in JK:
    if i < 4100:
        ListaJK.append(i)
for i in VJ:
    if i < 4100:
        ListaVJ.append(i)
        

for i in range(10):
    ListaJK.append(i)
    

MstJK = np.array(ListaJK)
MstVJ = np.array(ListaVJ)

len(MstVJ)


# In[7]:


fig, ax = plt.subplots()

ax.set(xlabel='Cor V-J', ylabel='Cor J-K')
ax.grid()

ax.scatter(Teff4, Teff7, s =1, c = 'black')
ax.scatter(MstVJ, MstJK, s =1, c = 'red')


plt.xlim(2000, 8000)
plt.ylim(2000, 8000)
plt.show()
plt.clf()

#fig.savefig("fig.png")


# In[8]:


JKarr = np.array(J_Ks)
VJarr = np.array(V_J)

JK = np.array(Teff7)
JKlst = JK.tolist()
VJ = np.array(Teff4)
VJlst = VJ.tolist()


ListaJK=[]
ListaVJ=[]


for i in JKlst:
    if i < 4100:
        ListaJK.append(JKlst.index(i))
for i in VJlst:
    if i < 4100:
        ListaVJ.append(VJlst.index(i))
        
LstmagJK = []
LstmagVJ = []
        
        
for i in ListaJK:
    x = JKarr[i]
    LstmagJK.append(x)
for i in ListaVJ:
    y = VJarr[i]
    LstmagVJ.append(y)    
    
for i in range(10):
    LstmagJK.append(i)
    


# In[9]:


fig, ax = plt.subplots()

ax.set(xlabel='Cor V-J', ylabel='Cor J-K')
ax.grid()

ax.scatter(V_J, J_Ks, s =1, c = 'black')
ax.scatter(LstmagVJ, LstmagJK, s =1, c = 'red')


plt.xlim(-1, 7)
plt.ylim(-0.25, 1.75)
plt.show()
plt.clf()

#fig.savefig("magfig.png")


# In[16]:


fig, ax = plt.subplots()

ax.set(xlabel='Teff NASA', ylabel='Teff Raul')
ax.grid()

#ax.scatter(t['st_teff'], Teff7, s=1)
ax.scatter(t['st_teff'], Teff1, s=1, c = 'red')
#ax.scatter(t['st_teff'], Teff5, s=1, c = 'green')
ax.scatter(t['st_teff'], Teff4, s=1, c = 'black')
plt.xlim(1500, 8000)
plt.ylim(1500, 8000)
plt.show()
plt.clf()


# In[ ]:




