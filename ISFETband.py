#! /usr/bin/python3

import matplotlib.pyplot as pylab
    
    
    
def isfetb_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_density, p_doping_density, n_doping_density, req_bend, T, pK1, pK2, steric):
    

    if steric==True:
        import ISFETsolverSteric as isfs
    else:
        import ISFETsolver as isfs
    
    tol = 1e-9
    err = 1.0
    print("Looking for fluid bias for band bending of",req_bend,"V at salt concentration",salt_conc,"M, pH",pH,"and surface site density",site_density,"1/cm2")
    
    def band_bend_eval(rbias):

        global space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC 
        
        space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC = isfs.isfet_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_density, p_doping_density, n_doping_density, rbias, T, pK1, pK2)

        band_bending = V_grid[0]-V_grid[end_o]
 
        return band_bending


    rbias = 0.0
    delta_b = 0.001
    while abs(err)>tol:
        bb2 = band_bend_eval(rbias+delta_b)
        bb1 = band_bend_eval(rbias)
        jac = (bb2-bb1)/delta_b
        rbias = rbias - (bb1-req_bend)*(1/jac)

        err = bb1-req_bend
        print("Current Band Bending=",err)

    print("Calculated fluid bias for",req_bend,"V bend is", rbias)
    print("Band Bending in Semiconductor at calculated potential =",bb1)

    return space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC, rbias 

#flatband_bias_finder(1e-3,1,1e14)

## Plot Potential
if False:
    #pylab.plot(space_grid[:-1]*1e6,D,label="D")
    pylab.plot(space_grid*1e6,V_grid,label="V")
    pylab.axvline(space_grid[end_s]*1e6,color='green')
    pylab.axvline(space_grid[end_o]*1e6,color='green')
    pylab.xlabel("Distance (um)")
    pylab.ylabel("Potential (V)")
    pylab.text(0.133, 0.05, 'Semiconductor', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.5, 0.05, 'Oxide', rotation='vertical',horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.833, 0.05, 'Electrolyte', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.legend()
    pylab.show() 

## Plot salt in electrolyte
if False:
    pylab.plot(space_grid*1e6,Na_conc*1e-3,label="I+")
    pylab.plot(space_grid*1e6,Cl_conc*1e-3,label="I-")
    pylab.axvline(space_grid[end_s]*1e6,color='green')
    pylab.axvline(space_grid[end_o]*1e6,color='green')
    pylab.xlabel("Distance (um)")
    pylab.ylabel("Concentration (M)")
    pylab.text(0.133, 0.05, 'Semiconductor', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.5, 0.05, 'Oxide', rotation='vertical',horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.833, 0.05, 'Electrolyte', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.legend()
    pylab.show()

## Plot protons and counterion concentration in electrolyte
if False:
    pylab.plot(space_grid*1e6,H_conc*1e-3,label="H+")
    pylab.plot(space_grid*1e6,OH_conc*1e-3,label="OH-")

    pylab.axvline(space_grid[end_s]*1e6,color='green')
    pylab.axvline(space_grid[end_o]*1e6,color='green')
    pylab.xlabel("Distance (um)")
    pylab.ylabel("Concentration (M)")
    pylab.text(0.133, 0.05, 'Semiconductor', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.5, 0.05, 'Oxide', rotation='vertical',horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.833, 0.05, 'Electrolyte', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.legend()
    pylab.show()

## Plot electrons and holes concentration in semiconductor
if False:
    pylab.plot(space_grid*1e6,h_conc*1e-6,label="h+")
    pylab.plot(space_grid*1e6,e_conc*1e-6,label="e-")
    pylab.axvline(space_grid[end_s]*1e6,color='green')
    pylab.axvline(space_grid[end_o]*1e6,color='green')
    pylab.xlabel("Distance (um)")
    pylab.ylabel("Concentration (cm^-3)")
    pylab.text(0.133, 0.05, 'Semiconductor', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.5, 0.05, 'Oxide',rotation='vertical', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.text(0.833, 0.05, 'Electrolyte', horizontalalignment='center',verticalalignment='center', transform=pylab.gca().transAxes)
    pylab.legend()
    pylab.show()
