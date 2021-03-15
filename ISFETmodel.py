import matrixGenerator as mg
import numpy as np
import matplotlib.pyplot as pylab
import shelve
import math
from datetime import datetime
import sys
import os
import progressbar



## Physical Constants
kb = 1.38064852e-23
el_c = 1.60217662e-19
perm_0 = 8.854178176e-12
Na = 6.022140857e23
Far = Na*el_c



def isfetmodel_1D(L_list, min_spacing,rel_perm_list, salt_conc_list, pH_list, surf_states_list, p_doping_density, n_doping_density, T=300.0, pK1=6.7, pK2=-1.9, mode=['N',0.0], steric=False, shelf_loc="./", title = "Untitled", force_new_calc=False, plots=[True]*8):

    import ISFETplotter as php
    filename = shelf_loc+"shelve_"+title+".out"
    if force_new_calc == False and os.path.isfile(filename):
        print("Shelf file Already exist. Not calculating again")
        php.isfetmodel_1D_plot(filename,title,plots)
        exit()

    if steric==True:
        import ISFETsolverSteric as isfet
    else:
        import ISFETsolver as isfet

    import ISFETband as isfetb
    import ISFETcharge as isfetc


    if shelf_loc[-1] != '/':
        shelf_loc = shelf_loc+"/"

    startTime = datetime.now()
    try:
        os.stat(shelf_loc)
    except:
        os.mkdir(shelf_loc)
      
    
    ## Converting python list to numpy arrays
    salt_conc_list = np.array(salt_conc_list)
    surf_states_list = np.array(surf_states_list)
    pH_list = np.array(pH_list)


    for k,salt_conc in enumerate(salt_conc_list):
        
        for a,site_den in enumerate(surf_states_list):
            
            for i,pH in enumerate(pH_list):  
                
                
                if mode[0] == 'N':
                    rbias = mode[1]
                    space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC = isfet.isfet_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_den, p_doping_density, n_doping_density, mode[1], T, pK1, pK2)
                
                elif mode[0] == 'B':
                    space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC, rbias = isfetb.isfetb_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_den, p_doping_density, n_doping_density, mode[1], T, pK1, pK2, steric)
                
                elif mode[0] == 'C':
                    space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC, rbias = isfetc.isfetc_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_den, p_doping_density, n_doping_density, mode[1], T, pK1, pK2, steric)
                
                
                
                
                if (k+a+i)==0:
                    ## Shelving Variables Initialization
                    tup1 = (salt_conc_list.size,surf_states_list.size,pH_list.size,space_grid[end_o:].size)
                    tup2 = (salt_conc_list.size,surf_states_list.size,pH_list.size)
                    tup3 = (salt_conc_list.size,surf_states_list.size,pH_list.size,space_grid[:end_s+1].size)
                    tup4 = (salt_conc_list.size,surf_states_list.size,pH_list.size,space_grid.size)

                    conc_posI_shelve = np.zeros(tup1)
                    conc_H_shelve = np.zeros(tup1) 
                    conc_h_shelve = np.zeros(tup3)
                    conc_negI_shelve = np.zeros(tup1) 
                    conc_OH_shelve = np.zeros(tup1) 
                    conc_e_shelve = np.zeros(tup3)
                    
                    V_grid_shelve = np.zeros(tup4) 

                    conc_SiO_shelve = np.zeros(tup2) 
                    conc_SiOH_shelve = np.zeros(tup2) 
                    conc_SiOH2_shelve = np.zeros(tup2)
                    srf_pot_shelve = np.zeros(tup2) 
                    rho_SC_shelve = np.zeros(tup2)
                    
                    rbias_shelve = np.zeros(tup2)
                
                ## Saving Results
                conc_posI_shelve[k][a][i] = conc_posI[end_o:]
                conc_H_shelve[k][a][i] = conc_H[end_o:]
                conc_h_shelve[k][a][i] = conc_h[:end_s+1]
                conc_negI_shelve[k][a][i] = conc_negI[end_o:]
                conc_OH_shelve[k][a][i] = conc_OH[end_o:]
                conc_e_shelve[k][a][i] = conc_e[:end_s+1]
                
                V_grid_shelve[k][a][i] = V_grid
                
                conc_SiO_shelve[k][a][i] = conc_SiO
                conc_SiOH_shelve[k][a][i] = conc_SiOH
                conc_SiOH2_shelve[k][a][i] = conc_SiOH2
                
                srf_pot_shelve[k][a][i] = V_grid[end_o]
                
                rho_SC_shelve[k][a][i] = rho_SC

                rbias_shelve[k][a][i] = rbias

                


                ## Status Printing
                print("Done with calculation for pH="+str(pH_list[i])+" SS =  "+str(surf_states_list[a])+" 1/cm^2"+" IC =  "+str(salt_conc_list[k])+" M. "+str((pH_list.size*salt_conc_list.size*surf_states_list.size)-(i+1)-(a*pH_list.size)-(k*pH_list.size*surf_states_list.size))+" calculations remaining")
                print("Surface Potential = "+str(V_grid[end_o]*1e3)+" mV \n")
            


    ## Saving Date File
    
    try:
        my_shelf = shelve.open(filename,'n')
    except:
        os.remove(filename)
        my_shelf = shelve.open(filename,'n')

    for key in ['conc_SiO_shelve','conc_SiOH_shelve','conc_SiOH2_shelve','srf_pot_shelve','conc_H_shelve','conc_posI_shelve','conc_h_shelve','conc_OH_shelve','conc_negI_shelve','conc_e_shelve','V_grid_shelve','pH_list','surf_states_list','salt_conc_list','space_grid','end_s','end_o','rho_SC_shelve','rbias_shelve']:
        try:
            my_shelf[key] = locals()[key]
        except TypeError:
            print('ERROR shelving: {0}'.format(key))
    my_shelf.close()
    

    ## Plotting Results
    php.isfetmodel_1D_plot(filename,title,plots)

    print("Total Execution Time :"+str(datetime.now() - startTime))



