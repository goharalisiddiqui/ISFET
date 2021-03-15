#! /usr/bin/python3

import matplotlib.pyplot as pylab
from datetime import datetime
import os
import sys
import shelve
import numpy as np
import progressbar
import matplotlib.ticker as ticker


def isfetmodel_1D_plot(shelf_file,title,plots):

    
    startTime = datetime.now()

    ## Some parameters
    Na = 6.022140857e23
    img_dpi = 400

    print("\n\n Plotting Now ...... \n")
    
    
    
    ### Plotting choices
    elec = plots[0]                    # Electrostatics on space grid
    SSforSS = plots[1]                # Three Species site density vs bulk pH for different total surface density
    SPforSS = plots[2]               # Surface potential vs bulk pH for different total surface density
    SpHforSS = plots[3]              # Surface pH vs bulk pH for different total surface density
    SRforSS = plots[4]              # Sensing ratio vs pH for different total surface density
    SRforIC = plots[5]            # Sensing ratio vs pH for different bulk ionic concentration
    SSforIC = plots[6]            # Three Species site density vs bulk pH for different bulk IC and highest total surface density
    SCforSS = plots[7]              # Total sheet charge of semiconductor vs bulk pH
    FBforSS = plots[8]


    out_folder = os.path.dirname(os.path.realpath(shelf_file)) + "/Results_"+title


    try:
        os.stat(out_folder)
    except:
        os.mkdir(out_folder) 


    ## Retrieving Shelved Variables

    filename = shelf_file
    
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        if SSforSS or SSforIC:
            conc_SiOH2 = np.array(my_shelf["conc_SiOH2_shelve"])
            conc_SiOH = np.array(my_shelf["conc_SiOH_shelve"])
            conc_SiO = np.array(my_shelf["conc_SiO_shelve"])
        if SPforSS or SRforSS or SRforIC:
            srf_pot = np.array(my_shelf["srf_pot_shelve"])
        if elec:
            space_grid = np.array(my_shelf["space_grid"])
            conc_posI = np.array(my_shelf["conc_posI_shelve"])
            conc_negI = np.array(my_shelf["conc_negI_shelve"])
            conc_H = np.array(my_shelf["conc_H_shelve"])
            conc_OH = np.array(my_shelf["conc_OH_shelve"])
        salt_conc_list = my_shelf["salt_conc_list"]
        V_grid = np.array(my_shelf["V_grid_shelve"])
        pH_list = my_shelf["pH_list"]
        surf_states_list = my_shelf["surf_states_list"]
        rho_SC = my_shelf["rho_SC_shelve"]
        rbias = my_shelf["rbias_shelve"]
        end_s = my_shelf["end_s"]
        end_o = my_shelf["end_o"]
    my_shelf.close()

    

    

    ## Progress Bar
    maxval = int(pH_list.size*surf_states_list.size*salt_conc_list.size)
    bar = progressbar.ProgressBar(maxval=maxval, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()




    ## pH values to plot
    pH_plot = np.unique(pH_list.astype(int))   # Which pH to plot the electrostatics for. Default plots for all the pH calculated


    ## Figues initialization for Sensing ration for different IC
    
    if SRforIC: 
        ICfig = []
        ICfigaxes = []
        for u in range(surf_states_list.size):
            fig, ax = pylab.subplots()
            ICfig.append(fig)
            ICfigaxes.append(ax)
    

    ## Figures for surface densities for the three species for differnet IC 
    if SSforIC:
        SSfig = []
        SSfigaxes = []
        for u in range(3):
            fig, ax = pylab.subplots()
            SSfig.append(fig)
            SSfigaxes.append(ax)

    
    for i in range(salt_conc_list.size):

        ## Making Directory for different salt concentrations results
        if elec or SSforSS or SPforSS or SpHforSS or SRforSS:
            ndir = out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]"
            try:
                os.stat(ndir)
            except:
                os.mkdir(ndir)


        ## Initializing plots
        if SSforSS:
            fig3, f3ax = pylab.subplots(nrows=3,ncols=1)
        if SPforSS:
            fig1, f1ax = pylab.subplots()
        if SpHforSS:
            fig4, f4ax = pylab.subplots()
        if SRforSS:
            fig5, f5ax = pylab.subplots()
        if SCforSS:
            figSCforSS, axSCforSS = pylab.subplots()
        if FBforSS:
            figFBforSS, axFBforSS = pylab.subplots()



        for j in range(surf_states_list.size):
            
            ###### Plotting Surface States
            if SSforSS:
                f3ax[0].plot(pH_list, conc_SiOH2[i][j], label="SS = "+str(surf_states_list[j])+" cm^-2")
                f3ax[1].plot(pH_list, conc_SiOH[i][j], label="SS = "+str(surf_states_list[j])+" cm^-2")
                f3ax[2].plot(pH_list, conc_SiO[i][j], label="SS = "+str(surf_states_list[j])+" cm^-2")
            

            ###### Plotting Surface densities of each species for differnet IC
            if j == 0 and SSforIC:
                
                SSfigaxes[0].plot(pH_list, conc_SiO[i][j], label="IC = "+str(salt_conc_list[i])+" M")
                SSfigaxes[1].plot(pH_list, conc_SiOH[i][j], label="IC = "+str(salt_conc_list[i])+" M")
                SSfigaxes[2].plot(pH_list, conc_SiOH2[i][j], label="IC = "+str(salt_conc_list[i])+" M")

            ###### Plotting Surface pH vs Bulk pH
            
            if SpHforSS:
                f4ax.plot(pH_list,(-np.log10(conc_H[i,j,:,0]*1e-3)), label="SS = "+str(surf_states_list[j])+" cm^-2")
            


            ###### Plotting Surface Potential Change vs pH for different SS
            if SRforSS:
                delta_srf_pot = (srf_pot[i][j][1:]-V_grid[i,j,1:,-1]-(srf_pot[i][j][:-1]-V_grid[i,j,:-1,-1]))/(pH_list[1:]-pH_list[:-1])
                f5ax.plot(pH_list[1:], delta_srf_pot*1e3, label="SS = "+str(surf_states_list[j])+" cm^-2")
            

            ###### Plotting Surface Potential Change vs pH for different IC
            if SRforIC:
                delta_srf_pot = (srf_pot[i][j][1:]-V_grid[i,j,1:,-1]-(srf_pot[i][j][:-1]-V_grid[i,j,:-1,-1]))/(pH_list[1:]-pH_list[:-1])
                ICfigaxes[j].plot(pH_list[1:], delta_srf_pot*1e3, label="IC = "+str(salt_conc_list[i])+" M")
            
            
            
            ###### Plotting Surface Potential vs pH
            if SPforSS:
                f1ax.plot(pH_list, (srf_pot[i][j]-V_grid[i,j,:,-1])*1e3, label="SS = "+str(surf_states_list[j])+" cm^-2")
            
            ###### Plotting Sheet Charge density vs pH
            if SCforSS:
                axSCforSS.plot(pH_list, rho_SC[i][j]*1e-4, label="SS = "+str(surf_states_list[j])+" cm^-2")
            
            
            ###### Plotting Fluid bias vs pH
            if FBforSS:
                axFBforSS.plot(pH_list, rbias[i][j]*1e3, label="SS = "+str(surf_states_list[j])+" cm^-2")
            

            ##### Making directories to store electrostatics plots
            if elec:
                ndir = out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]"+"/"+"SC_"+str(surf_states_list[j])+"[cm^-2]"
                try:
                    os.stat(ndir)
                except:
                    os.mkdir(ndir)
            
            
            
            for l in range(pH_list.size):
                
                if (pH_list[l] in pH_plot) and elec:

                    ## Correct space grid of fluid only
                    x = space_grid[end_o:]-space_grid[end_o]

                    ##### Searching for a suitable length to plot the Electrostatics
                    zeropoint = x.size-1
                    for ind in range(x.size):
                        if abs(V_grid[i][j][l][end_o+ind]-V_grid[i][j][l][-1]) < 1e-12:
                            zeropoint = ind
                            break
                    

                    ## Plots Initialization
                    fig2, f2ax = pylab.subplots(nrows=3, ncols=1)
                    fig2.tight_layout()


                    ###### Plotting Electrostatics
                    f2ax[0].plot(x[:zeropoint+1]*1e9, V_grid[i][j][l][end_o:end_o+zeropoint+1]*1e3, label="V")
                    f2ax[0].annotate("Bias = %.2f mV"%(V_grid[i][j][l][end_o+zeropoint]*1e3),(x[zeropoint]*1e9,V_grid[i][j][l][end_o+zeropoint]*1e3),(-300,-200),textcoords='offset pixels',fontsize='small')

                    f2ax[1].semilogy(x[:zeropoint+1]*1e9, conc_posI[i][j][l][:zeropoint+1]*1e-3, label="Positive salt ion")
                    f2ax[1].semilogy(x[:zeropoint+1]*1e9, conc_negI[i][j][l][:zeropoint+1]*1e-3, label="Negative salt ion")

                    
                    f2ax[2].plot(x[:zeropoint+1]*1e9, -np.log10(conc_H[i][j][l][:zeropoint+1]*1e-3), label="pH")
                    f2ax[2].plot(x[:zeropoint+1]*1e9, -np.log10(conc_OH[i][j][l][:zeropoint+1]*1e-3), label="pOH")

                    ##### Setting plot attributes for Electrostatics
                    f2ax[0].set_title("ISFET Electrostatics plot for pH = "+str(pH_list[l])+", SC = "+str(surf_states_list[j])+" cm^-2  and IC = "+str(salt_conc_list[i])+" M ",fontsize='small')
                    f2ax[0].set_ylabel("Potential (mV)",fontsize='small')
                    f2ax[0].legend()
                    
                    f2ax[1].set_ylabel("Concentration (M)",fontsize='small')
                    f2ax[1].legend()
                    
                    f2ax[2].set_xlabel("Distance from oxide surface (nm)")
                    f2ax[2].set_ylabel("pH/pOH",fontsize='small')
                    #f2ax[2].yaxis.set_major_locator(ticker.MultipleLocator(1.0))
                    f2ax[2].legend()
                    

                    ##### Saving and closing Electrostatics plots
                    fig2.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"SC_"+str(surf_states_list[j])+"[cm^-2]/Electrostatics_pH_"+str(pH_list[l])+".jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
                    pylab.close(fig2)
                
                bar.update((l+1)+(j*len(pH_list))+(i*len(pH_list)*len(surf_states_list)))


        ##### Setting plot attributes and saving Surface potential plots
        if SPforSS:
            f1ax.set_xlabel("pH",fontsize=14)
            f1ax.set_ylabel("Surface Potential difference with fluid bias (mV)",fontsize=14)
            f1ax.set_title("Surface PD with fluid bias for IC = "+str(salt_conc_list[i])+" M ")
            f1ax.set_xticks(pH_plot)
            #f1ax.set_ylim(-300, 10)
            f1ax.tick_params(labelsize=14)
            f1ax.minorticks_on()
            f1ax.grid(True)
            f1ax.legend()
            fig1.tight_layout()
            fig1.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Surface_Potential_difference.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig1)
        
        ##### Setting plot attributes and saving Sheet charge plots
        if SCforSS:
            axSCforSS.set_xlabel("pH",fontsize=14)
            axSCforSS.set_ylabel("Sheet charge density(C/cm2)",fontsize=14)
            axSCforSS.set_title("Sheet charge density for IC = "+str(salt_conc_list[i])+" M ")
            axSCforSS.set_xticks(pH_plot)
            #axSCforSS.set_ylim(-300, 10)
            axSCforSS.tick_params(labelsize=14)
            axSCforSS.minorticks_on()
            axSCforSS.grid(True)
            axSCforSS.legend()
            figSCforSS.tight_layout()
            figSCforSS.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Sheet_Charge.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(figSCforSS)
        
        
        
        
        
        #### Setting attributes and saving Surface States plots 
        if SSforSS:
            f3ax[0].set_title("Surface States for IC = "+str(salt_conc_list[i])+" M ")
            f3ax[0].set_title("SC of SiOH2+ (cm^-2)",fontsize=14)
            f3ax[1].set_title("SC of SiOH (cm^-2)",fontsize=14)
            f3ax[2].set_title("SC of SiO- (cm^-2)",fontsize=14)
            f3ax[2].set_xlabel("pH")
            for axx in f3ax:
                axx.set_xticks(pH_plot)
                axx.grid(True)
                axx.tick_params(labelsize=14)
                #axx.set_ylim(0, 1e14)
                axx.minorticks_on()
                axx.legend(prop={'size': 6})
            fig3.tight_layout()
            fig3.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Surface_States.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig3)
        
       


        #### Setting attributes and saving Surface pH plots
        if SpHforSS:
            f4ax.plot(pH_list,pH_list, label="Reference line x=y", linestyle='dashed')
            f4ax.set_xlabel("Bulk pH",fontsize=14)
            f4ax.set_ylabel("Surface_pH",fontsize=14)
            f4ax.set_title("Surface pH for IC = "+str(salt_conc_list[i])+" M ")
            f4ax.set_xticks(pH_plot)
            f4ax.set_yticks(pH_plot)
            f4ax.tick_params(labelsize=14)
            f4ax.minorticks_on()
            f4ax.grid(True)
            f4ax.legend()
            fig4.tight_layout()
            fig4.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Surface_pH.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig4)



        #### Setting attributes and saving the Sensing ratio plots
        if SRforSS:
            f5ax.set_xlabel("pH",fontsize=14)
            f5ax.set_ylabel("Sensing Ratio (mV/pH)",fontsize=14)
            f5ax.set_title("Sensing ratio for IC = "+str(salt_conc_list[i])+" M ")
            f5ax.set_xticks(pH_plot)
            #f5ax.set_ylim(-40, 1)
            f5ax.tick_params(labelsize=14)
            f5ax.minorticks_on()
            f5ax.grid(True)
            f5ax.legend()
            fig5.tight_layout()
            fig5.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Sensing_ratio.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig5)
        

        ##### Setting plot attributes and saving Fluid bias plots
        if FBforSS:
            axFBforSS.set_xlabel("pH",fontsize=14)
            axFBforSS.set_ylabel("Fluid bias  (mV)",fontsize=14)
            axFBforSS.set_title("Fluid bias for IC = "+str(salt_conc_list[i])+" M ")
            axFBforSS.set_xticks(pH_plot)
            #axFBforSS.set_ylim(-300, 10)
            axFBforSS.tick_params(labelsize=14)
            axFBforSS.minorticks_on()
            axFBforSS.grid(True)
            axFBforSS.legend()
            figFBforSS.tight_layout()
            figFBforSS.savefig(out_folder+"/"+"IC_"+str(salt_conc_list[i])+"[M]/"+"Fluid_Bias.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(figFBforSS)
        
    #### Setting attributes and saving the Sensing ration for different IC plots
    if SRforIC:
        for figax in ICfigaxes:
            #### Plot attributes for the sensing change with IC plot
            figax.set_xlabel("pH",fontsize=14)
            figax.set_ylabel("Sensing ratio (mV/pH)",fontsize=14)
            figax.set_title("Sensing ratio for different IC")
            figax.set_xticks(pH_plot)
            figax.tick_params(labelsize=14)
            figax.minorticks_on()
            figax.grid(True)
            #figax.set_ylim(-40, 1)
            figax.legend()

        for i,fig in enumerate(ICfig):
            ##### Saving and Closing the last plot
            fig.tight_layout()
            fig.savefig(out_folder+"/"+"SC_"+str(surf_states_list[i])+"cm^-2"+"_Sensing_ratio.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig)
   




    #### Setting attributes and saving the Surface density for different IC plots
    if SSforIC:
        for ind,figax in enumerate(SSfigaxes):
            #### Plot attributes for the sensing change with IC plot
            figax.set_xlabel("pH",fontsize=14) 
            figax.set_title("Surface states density for different IC")
            if ind == 0:
                figax.set_ylabel("[SiO-] in cm^-2",fontsize=14)
            elif ind == 1:
                figax.set_ylabel("[SiOH] in cm^-2",fontsize=14)
            elif ind == 2:
                figax.set_ylabel("[SiOH2+] in cm^-2",fontsize=14)
            figax.set_xticks(pH_plot)
            figax.tick_params(labelsize=14)
            figax.minorticks_on()
            figax.grid(True)
            #figax.set_ylim(-40, 1)
            figax.legend()

        for i,fig in enumerate(SSfig):
            ##### Saving and Closing the last plot
            fig.tight_layout()
            if i == 0:
                fig.savefig(out_folder+"/"+"SiO_Surface_density.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            elif i == 1:
                fig.savefig(out_folder+"/"+"SiOH_Surface_density.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            elif i == 2:
                fig.savefig(out_folder+"/"+"SiOH2_Surface_density.jpg",dpi=img_dpi,bbox_inches="tight",quality=95)
            pylab.close(fig)
    
    
    
    
    bar.finish()
    print("Plotting Time :"+str(datetime.now() - startTime))
