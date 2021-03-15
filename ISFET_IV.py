#! /usr/bin/python3

import matplotlib.pyplot as pylab
from datetime import datetime
import os
import sys
import shelve
import numpy as np
import progressbar
import matplotlib.ticker as ticker
import glob
import time

    

## Some parameters
Na = 6.022140857e23
img_dpi = 400
upd = True




pH_plot = [1.0,5.0,9.0]#0.1*np.array(range(10,100,10)) #np.unique(pH_list.astype(int))   # Which pH to plot the electrostatics for. Default plots for all the pH calculated

figFB, axFB = pylab.subplots()

## Retrieving Shelved Variables
file_list = sorted(glob.glob('./Results_ISFET/shelve_ISFET_IV_0*.out'))
for ind,filename in enumerate(file_list):
    if upd:
        rbias = [[] for x in range(len(file_list))]
        rho_SC = [[] for x in range(len(file_list))]
        upd = False
    my_shelf = shelve.open(filename)
    for i,p in enumerate(my_shelf["pH_list"]):
        if p in pH_plot:
            rho_SC[ind].append(my_shelf["rho_SC_shelve"][0][0][i])
            rbias[ind].append(my_shelf["rbias_shelve"][0][0][i])
    my_shelf.close()

rho_SC = np.transpose(np.array(rho_SC))
rbias = np.transpose(np.array(rbias))


for ind,pH in enumerate(pH_plot):
    axFB.plot(rbias[ind], rho_SC[ind]*1e2,'o-', label="pH = %d"%(pH))

    

    
#axFB.set_yscale('log')
##### Setting plot attributes and saving Fluid bias plots
axFB.set_xlabel("Fluid bias (V)",fontsize=14)
axFB.set_ylabel("Sheet Charge density (uC/cm2)",fontsize=14)
axFB.set_title("Transfer Characteristic for different pH values")
axFB.tick_params(labelsize=16)
axFB.minorticks_on()
axFB.grid(True)
axFB.legend(fontsize=14)
figFB.tight_layout()
pylab.show() 
#pylab.savefig("IV.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
#pylab.close(figFB)
        
