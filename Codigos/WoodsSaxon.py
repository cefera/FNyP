#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Física Nuclear y de Partículas

@author: cesar fernandez ramirez

version: junio 2024
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

#%%

def radio(r0,A):
    return  r0*(A**(1/3))

def densidadnuclear(r0,A,a,r):
    R = radio(r0,A)
    return 1./(1.+np.exp((r-R)/a))

def WoodsSaxon(r0,A,a,r):
    R = radio(r0,A)
    return -1./(1.+np.exp((r-R)/a))


fig = plt.figure(figsize=(5,5))
r = np.arange(0.,10.,0.1)

plt.ylabel(r'$\rho/\rho_0$',fontsize=15)
plt.xlabel(r'$r$ (fm)',fontsize=15)
plt.xlim((0.,10.))
plt.ylim((0.,1.))
plt.plot(r,densidadnuclear(1.2,12,0.55,r),'-',c=jpac_color[0],label='A=12')
plt.plot(r,densidadnuclear(1.2,40,0.55,r),'-',c=jpac_color[1],label='A=40')
plt.plot(r,densidadnuclear(1.2,208,0.55,r),'-',c=jpac_color[2],label='A=208')
plt.vlines(radio(1.2,12),0.,1.,colors=jpac_color[0], lw=1., linestyles='dashed')
plt.vlines(radio(1.2,40),0.,1.,colors=jpac_color[1], lw=1., linestyles='dashed')
plt.vlines(radio(1.2,208),0.,1.,colors=jpac_color[2], lw=1., linestyles='dashed')
plt.text(radio(1.2,12)+0.1,0.1,r'$R$',c=jpac_color[0],fontsize=15)
plt.text(radio(1.2,40)+0.1,0.1,r'$R$',c=jpac_color[1],fontsize=15)
plt.text(radio(1.2,208)+0.1,0.1,r'$R$',c=jpac_color[2],fontsize=15)
plt.legend(loc='upper right',ncol=1,frameon=True)
plt.show()    
fig.savefig("densidadnuclear.pdf", bbox_inches='tight')


fig = plt.figure(figsize=(5,5))
r = np.arange(0.,10.,0.1)
plt.ylabel(r'$V/V_0$',fontsize=15)
plt.xlabel(r'$r$ (fm)',fontsize=15)
plt.xlim((0.,10.))
plt.ylim((-1.,0.2))
plt.plot(r,WoodsSaxon(1.2,40,0.55,r),'-',c=jpac_color[0])
plt.hlines(-0.90,0.,10.,colors=jpac_color[1], lw=1., linestyles='dashed',label='Fundamental')
plt.hlines(-0.80,0.,10.,colors=jpac_color[2], lw=1., linestyles='dashed',label='Primer excitado')
plt.hlines(-0.75,0.,10.,colors=jpac_color[3], lw=1., linestyles='dashed',label='Segundo excitado')
plt.hlines(0.,0.,10.,colors=jpac_color[9], lw=1., linestyles='solid')
plt.text(0.5,-0.1,r'Discreto, estados ligados',c=jpac_color[0],fontsize=10)
plt.text(0.5,0.1,r'Continuo',c=jpac_color[1],fontsize=10)
plt.legend(loc='center right',ncol=1,frameon=True)
plt.show()    
fig.savefig("WoodsSaxon.pdf", bbox_inches='tight')

#%%

def coulomb(Z,r):
    return -2*Z/r

def barrera(l,r):
    return l*(l+1)/(r*r)

def potencial(Z,l,r):
    return coulomb(Z,r) + barrera(l,r)

def En(ry,nr,l):
    n = nr + l + 1
    return -ry/(n**2)

a0,ry  = 5.291772105, 13.6058
Z, lmax, nmax =1, 2, 5

fig = plt.figure(figsize=(5,5))
rmin, rmax, rstep = 0.1, 15, 0.01
r = np.arange(rmin,rmax,rstep)
plt.xlim((rmin,rmax))
plt.ylim((-20.,5))
plt.ylabel(r'$V$ (eV)',fontsize=15)
plt.xlabel(r'$r$ ($a_0$)',fontsize=15)

plt.hlines(0,rmin,rmax,colors=jpac_color[9], lw=1., linestyles='solid', alpha=0.5)
for l in range(lmax+1):
    texto = '$\ell=$' +  str(l)
    plt.plot(r,potencial(Z,l,r)*ry,'-',c=jpac_color[l],label=texto)
    plt.hlines(En(ry,0,l),rmin,rmax,colors=jpac_color[l], lw=1., linestyles='dashed')
    
    for n in range(1,nmax+1):
        plt.hlines(En(ry,n,l),rmin,rmax,colors=jpac_color[l], lw=1., linestyles='dotted', alpha=0.5)

plt.legend(loc='lower right',ncol=1,frameon=True)
plt.show()
fig.savefig("hidrogeno.pdf", bbox_inches='tight')






