#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Física Nuclear y de Partículas

@author: cesar fernandez ramirez

version: junio 2024


Based on Christian Hill's https://scipython.com/blog/nuclear-binding-energies-1/ 
"""

##############################################################################   
#
#   Libraries
#
##############################################################################   

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###########################################################
#   JPAC color style
###########################################################

jpac_blue   = "#1F77B4"; jpac_red    = "#D61D28";
jpac_green  = "#2CA02C"; jpac_orange = "#FF7F0E";
jpac_purple = "#9467BD"; jpac_brown  = "#8C564B";
jpac_pink   = "#E377C2"; jpac_gold   = "#BCBD22";
jpac_aqua   = "#17BECF"; jpac_grey   = "#7F7F7F";

jpac_color = [jpac_blue, jpac_red, jpac_green, 
              jpac_orange, jpac_purple, jpac_brown,
              jpac_pink, jpac_gold, jpac_aqua, jpac_grey, 'black' ];

jpac_axes = jpac_color[10]
dashes = 40*'-'

def SEMF(Z,N):
    aV, aS, aC, aA, delta = 15.75, 17.8, 0.711, 23.7, 11.18
    Z, N = np.atleast_1d(Z), np.atleast_1d(N)
    A = Z + N
    sgn = np.zeros(Z.shape)
    sgn[(Z%2) & (N%2)] = -1.
    sgn[~(Z%2) & ~(N%2)] = +1.
    return aV - aS / A**(1./3.) - aC * Z**2. / A**(4./3.) - aA * (A-2.*Z)**2./A**2. + sgn * delta/A**(3./2.)

file = 'mass_1.mas20.txt'

df = pd.read_fwf(file, usecols=(2,3,4,11),
              widths=(1,3,5,5,5,1,3,4,1,13,11,11,9,1,2,11,9,1,3,1,12,11,1),
              skiprows=39, header=None,
              index_col=False)
df.columns = ('N', 'Z', 'A', 'avEbind')
df['avEbind'] = pd.to_numeric(df['avEbind'], errors='coerce')
df = df.dropna()
# Also convert from keV to MeV.
df['avEbind'] /= 1000
gdf = df.groupby('A')
maxavEbind = gdf.apply(lambda t: t[t.avEbind==t.avEbind.max()])
maxavEbind['Eapprox'] = SEMF(maxavEbind['Z'], maxavEbind['N'])

fig = plt.figure(figsize=(8,5))
plt.ylabel(r'$E_B/A$ (MeV)',fontsize=15)
plt.xlabel(r'$A$',fontsize=15)
plt.ylim((1.,9.))
plt.xlim((0.,280.))
plt.hlines(8.09670568181818,0.,280.,colors=jpac_color[9], lw=1., linestyles='dashed')
plt.hlines(7.25,0.,280.,colors=jpac_color[9], lw=1., linestyles='dashed')
plt.vlines(62,0.,maxavEbind['avEbind'].max(),colors=jpac_color[3], lw=1., linestyles='dashed')
plt.scatter(maxavEbind['A'], maxavEbind['avEbind'],marker='o',s=5,c=jpac_color[0],label='AME2020')
plt.plot(maxavEbind['A'], maxavEbind['avEbind'],'-',lw=1,c=jpac_color[0])
plt.text(215,8.2,r'$\bar{E}_B=8,1$ MeV',c=jpac_color[9],fontsize=15)
plt.text(64,2,r'$E_B(^{62}_{28}$Ni$)=8,7945$ MeV',c=jpac_color[9],fontsize=15)
plt.text(210,6.6,r'$E_B\simeq 7,25$ MeV',c=jpac_color[9],fontsize=15)
plt.text(210,6.1,r'$A=270$',c=jpac_color[9],fontsize=15)
plt.legend(loc='center right',ncol=1,frameon=True)
plt.show()    
fig.savefig("BindingEnergy.pdf", bbox_inches='tight')

fig = plt.figure(figsize=(8,5))
plt.ylabel(r'$E_B/A$ (MeV)',fontsize=15)
plt.xlabel(r'$A$',fontsize=15)
plt.ylim((1.,9.))
plt.scatter(maxavEbind['A'], maxavEbind['avEbind'],marker='o',s=5,c=jpac_color[0],label='AME2020')
plt.plot(maxavEbind['A'], maxavEbind['avEbind'],'-',lw=1,c=jpac_color[0])
plt.plot(maxavEbind['A'], maxavEbind['Eapprox'],'-',lw=1,c=jpac_color[1],label='Bethe-Weizsacker')
plt.legend(loc='center right',ncol=1,frameon=True)
plt.show()    
fig.savefig("BindingEnergy2.pdf", bbox_inches='tight')

fig = plt.figure(figsize=(8,5))
plt.ylabel(r'$E_B/A-E_B(Z,A)/A$ (MeV)',fontsize=15)
plt.xlabel(r'$A$',fontsize=15)
plt.ylim((-1.,1.))
plt.xlim((7.,280.))
plt.scatter(maxavEbind['A'], maxavEbind['avEbind']-maxavEbind['Eapprox'],marker='o',s=5,c=jpac_color[0],label='Exp-Teo')
plt.plot(maxavEbind['A'], maxavEbind['avEbind']-maxavEbind['Eapprox'],'-',lw=1,c=jpac_color[0])
plt.legend(loc='upper right',ncol=1,frameon=True)
plt.show()    
fig.savefig("BindingEnergy3.pdf", bbox_inches='tight')

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(projection='3d')
ax.set_xlabel('Z')
ax.set_ylabel('N')
ax.set_zlabel(r'$E_B$ (MeV)')
ax.scatter(df['Z'],df['N'],df['avEbind'],marker='o',s=1,c=jpac_color[0])
plt.show()    
fig.savefig("BindingEnergy3D.pdf", bbox_inches='tight')


print(' '); print(dashes); print(' ')
print('Average binding energy of all nuclei:',df['avEbind'].mean())
print('Standard deviation:',df['avEbind'].std())
print(' '); print(dashes); print(' ')
print('Nucleus with the highest binding energy: Exp')
print(maxavEbind[ maxavEbind['avEbind'] == maxavEbind['avEbind'].max()])
print(' ')
print('Average binding energy:',maxavEbind['avEbind'].mean())
print(' '); print(dashes); print(' ')
print('Nucleus with the highest binding energy: Theo')
print(maxavEbind[ maxavEbind['Eapprox'] == maxavEbind['Eapprox'].max()])
