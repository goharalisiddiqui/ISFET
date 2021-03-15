#! /usr/bin/python3

import matplotlib.pyplot as pylab
    
    
    
def isfetc_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_density, p_doping_density, n_doping_density, req_char, T, pK1, pK2, steric):
    
    
    
    if steric==True:
        import ISFETsolverSteric as isfs
    else:
        import ISFETsolver as isfs
    
    
    tol = 1e-9
    err = 1.0
    print("Looking for fluid bias for sheet charge density of", req_char,"C/cm2 at salt concentration",salt_conc,"M, pH",pH,"and surface site density",site_density,"1/cm2")
    req_charge = req_char*1e4
    def sheet_charge_eval(rbias):

        global space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC 
        
        space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC = isfs.isfet_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_density, p_doping_density, n_doping_density, rbias, T, pK1, pK2)

        return rho_SC


    rbias = 0.0
    delta_b = 0.001
    while abs(err)>tol:
        bb2 = sheet_charge_eval(rbias+delta_b)
        bb1 = sheet_charge_eval(rbias)
        jac = (bb2-bb1)/delta_b
        y = (bb1-req_charge)*(1.0/jac)
        y = min(max(y,-0.2),0.2)
        rbias = rbias - y 

        err = bb1-req_charge
        print("Current Sheet Charge difference=",err)

    print("Calculated fluid bias for",req_char,"C/cm2 sheet charge is", rbias)
    print("Sheet charge in Semiconductor at calculated potential =",bb1*1e-4)

    return space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO, conc_SiOH, conc_SiOH2, rho_SC, rbias 

