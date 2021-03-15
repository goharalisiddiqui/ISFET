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


    

## Some parameters
Na = 6.022140857e23
img_dpi = 400






figFB, axFB = pylab.subplots()

## Retrieving Shelved Variables

for filename in sorted(glob.glob('./Results_ISFET/shelve_ISFET_FB_0*_stericSC.out')):
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        pH_list = my_shelf["pH_list"]
        rho_SC = my_shelf["rho_SC_shelve"]
        rbias = my_shelf["rbias_shelve"]
    my_shelf.close()
    ###### Plotting Fluid bias vs pH
    axFB.plot(pH_list, (rbias[0][0]-rbias[0][0][-1])*1e3, label="SC = %.2E C/cm^-2"%(rho_SC[0][0][0]*1e-4 if abs(rho_SC[0][0][0]*1e-4)> 1e-13 else 0.0))

    

pH_plot = np.unique(pH_list.astype(int))   # Which pH to plot the electrostatics for. Default plots for all the pH calculated
    

        ##### Setting plot attributes and saving Fluid bias plots
axFB.set_xlabel("pH",fontsize=14)
axFB.set_ylabel("Fluid bias  (mV)",fontsize=14)
axFB.set_title("Fluid bias for differnet sheet charge densities")
axFB.set_xticks(pH_plot)
axFB.tick_params(labelsize=14)
axFB.minorticks_on()
axFB.grid(True)
axFB.legend(fontsize=8)
figFB.tight_layout()
figFB.savefig("Fluid_Bias.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
pylab.close(figFB)
        
