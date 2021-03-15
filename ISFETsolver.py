import numpy as np
import matrixGenerator as mg
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve as spsolve
from scipy.linalg import solve
from scipy.linalg  import norm
import math
import log_grid_gen as lgg






def isfet_1D(L_list, min_spacing, rel_perm_list, salt_conc, pH, site_density, p_doping_density, n_doping_density, rbias, T, pK1, pK2):
    """
    Solve 1D ISFET model with supplied parameters and return the space grid, potential and all the concentrations(see return statement)

    Arguments:
    L_list = A list contaning the length in meters of the semiconductor, oxide and fluid region in that order
    min-spacing = The minimum spacing in meters for the space grid
    rel_perm_list = A list of relative permittivities of the semiconductor, oxide and fluid region in that order
    salt_conc = The salt concentration in Molar
    pH = The pH of the fluid
    site_density = The surface sites density of the oxide surface in cm^-2
    p_doping_density = Concentration in cm^-3 of acceptor dopants in the semiconductor
    n_doping_density = Concentration in cm^-3 of donor dopants in the semiconductor
    rbias = Voltage bias in Volts applied to the fluid
    T = Absolute temperature
    K1 = Equilibrium constant for trap state SiOH [m^3/moles.sec]
    K2 = Equilibrium constant for trap state SiOH2+ [m^3/moles.sec] 

    """

    ## Physical Constants
    kb = 1.38064852e-23
    el_c = 1.60217662e-19
    perm_0 = 8.854178176e-12
    Na = 6.022140857e23
    Far = Na*el_c


    ## Computational Parameters
    tol = 1e-8

    ## Physical Parameters
    K1 = math.pow(10,pK1-3)
    K2 = math.pow(10,pK2-3)
    perm_s = rel_perm_list[0]*perm_0                                                
    perm_o = rel_perm_list[1]*perm_0
    perm_f = rel_perm_list[2]*perm_0
    bulk_salt = salt_conc*1e3
    bulk_proton = math.pow(10,3-pH)
    site_den = (1e4/Na)*site_density
    doping_n = (n_doping_density)*1e6
    doping_p = (p_doping_density)*1e6
    int_den = (9.65e9)*1e6

    ## Derived Parameters
    V_T = (kb*T)/el_c


    ## Applying bias to the fluid
    fermi_fluid_pos = math.exp(rbias/V_T)
    fermi_fluid_neg = math.exp(-rbias/V_T)




    ### Initializations

    ## Space Grid
    min_spacing = 1e-10 #Minimum spacing for the space grids 
    L_s = L_list[0]
    L_o = L_list[1]
    L_f = L_list[2]

    space_grid_s = L_s-np.flip(lgg.loggrid(L_s,min_spacing),0)

    NPSD_o = int(L_o/min_spacing)+1
    space_grid_o = np.linspace(L_s,L_s+L_o,NPSD_o)[1:]

    space_grid_f = (L_s+L_o+lgg.loggrid(L_f,min_spacing))[1:]


    space_grid = np.append(np.append(space_grid_s,space_grid_o),space_grid_f)

    h = space_grid[1:] - space_grid[:-1]
    h_prev = np.append([0.0],h[:-1])

    L = L_s+L_o+L_f
    NPSD_s = space_grid_s.size
    NPSD_o = space_grid_o.size
    NPSD_f = space_grid_f.size
    NPSD = space_grid.size


    ## Permitivity Grid
    perm_grid = np.append(np.append(np.full(NPSD_s-1,perm_s),np.full(NPSD_o,perm_o)),np.full(NPSD_f+1,perm_f))
    perm_grid[-1] = 0.0 


    ## Matrix setup
    cap_mat = mg.get_capMat_1D_cartesian(space_grid, perm_grid,'N','D')

    ## Final indexes of regions
    end_s = NPSD_s-1
    end_o = NPSD_s+NPSD_o-1


    ##Potential grid
    V_grid = np.zeros(NPSD)


    ## Defect density
    conc_neg_D = np.zeros(NPSD)
    conc_pos_D = np.zeros(NPSD)         
    conc_neg_D[:end_s+1].fill(doping_p)
    conc_pos_D[:end_s+1].fill(doping_n)

    ## Trap states
    conc_SiO = np.zeros(NPSD)
    conc_SiOH = np.zeros(NPSD)
    conc_SiOH2 = np.zeros(NPSD)

    ## Charged species
    conc_posI = np.zeros(NPSD)
    conc_H = np.zeros(NPSD)
    conc_e = np.zeros(NPSD)
    conc_negI = np.zeros(NPSD)
    conc_OH = np.zeros(NPSD)
    conc_h = np.zeros(NPSD)
    
    
    delta = 1.0
    num_it = 0
    while delta > tol:


        ## Salt ions Charges
        conc_posI[end_o:] = bulk_salt*fermi_fluid_pos*np.exp(-V_grid[end_o:]/V_T)
        conc_negI[end_o:] = bulk_salt*fermi_fluid_neg*np.exp(V_grid[end_o:]/V_T)

        ## Protons and counter charge
        conc_H[end_o:] = bulk_proton*fermi_fluid_pos*np.exp(-V_grid[end_o:]/V_T)
        conc_OH[end_o:] = bulk_proton*fermi_fluid_neg*np.exp(V_grid[end_o:]/V_T)
       

        ## Electrons and holes
        conc_h[:end_s+1] = int_den*np.exp(-V_grid[:end_s+1]/V_T) 
        conc_e[:end_s+1] = int_den*np.exp(V_grid[:end_s+1]/V_T)

        
        ## Calculating surface states for current protons concentration at the surface
        M = np.array([[K1*conc_H[end_o],-1.0,0],[0,K2*conc_H[end_o],-1.0],[1,1,1]])
        C = np.array([0,0,site_den])
        X = solve(M,C)
        conc_SiO[end_o] = X[0]
        conc_SiOH[end_o] = X[1]
        conc_SiOH2[end_o] = X[2]

        ## Total charge
        rho_total = (Far*(conc_H + conc_posI - conc_OH - conc_negI - (conc_SiO/h[end_o]) + (conc_SiOH2/h[end_o]))) + (el_c*(conc_h - conc_e + conc_pos_D - conc_neg_D)) 
        
        ## Boundary conditions for charge density   
        rho_total[-1] = rbias                                                     #Dirichlet Condition for right boundary

        ## Total Function
        F_tot = (cap_mat.dot(V_grid)) - rho_total

        ## Calculating derivative of trap states for Jacobian
        M2 = np.array([[K1*conc_H[end_o],-1.0,0],[0,K2*conc_H[end_o],-1.0],[1,1,1]])
        C2 = np.array([-(K1*(-1.0/V_T)*conc_H[end_o]*conc_SiO[end_o]),-(K2*(-1.0/V_T)*conc_H[end_o]*conc_SiOH[end_o]),0])
        X2 = solve(M2,C2)
        Jac_SiO = X2[0]
        Jac_SiOH = X2[1]
        Jac_SiOH2 = X2[2]
        
        
        ## Jacobian
        P = (Far*(((-1.0/V_T)*(conc_H+conc_posI)) - ((1.0/V_T)*(conc_OH+conc_negI))))  +  (el_c*(((-1.0/V_T)*conc_h) - ((1.0/V_T)*conc_e)))
        P[end_o] = P[end_o] - ((Far*Jac_SiO)/h[end_o]) + ((Far*Jac_SiOH2)/h[end_o])
        P[-1] = 0.0 # Dirichtel boundary at right end
        
        Jac = cap_mat - (sparse.dia_matrix((P, [0]), shape=(NPSD, NPSD)))
        


        ## Solution
        y = spsolve(Jac, F_tot) ## Current Correction
        

        ## Limiting the step size in the first 50 iterations to avoid large fluctuations
        for m,j in enumerate(y):
            if j>V_T and num_it<50:
                y[m] = V_T
            elif j<-V_T and num_it<50:
                y[m] = -V_T

        
        ## Update potential
        V_grid = V_grid - y


        ## Calculate norm
        delta = norm(y)
        #print("Current Residual =", delta)


        ## Update iteration counter
        num_it += 1


    ### Calculating concentrations for final potential

    ## Salt ions Charges
    conc_posI[end_o:] = bulk_salt*fermi_fluid_pos*np.exp(-V_grid[end_o:]/V_T)
    conc_negI[end_o:] = bulk_salt*fermi_fluid_neg*np.exp(V_grid[end_o:]/V_T)

    ## Protons and counter charge
    conc_H[end_o:] = bulk_proton*fermi_fluid_pos*np.exp(-V_grid[end_o:]/V_T)
    conc_OH[end_o:] = bulk_proton*fermi_fluid_neg*np.exp(V_grid[end_o:]/V_T)
   

    ## Electrons and holes
    conc_h[:end_s+1] = int_den*np.exp(-V_grid[:end_s+1]/V_T) 
    conc_e[:end_s+1] = int_den*np.exp(V_grid[:end_s+1]/V_T)


    ## Calculating Sheet charge density
    rho_SC = 0.0
    for ind,point in enumerate(space_grid[:end_s+1]):
        rho_SC += rho_total[ind]*((h[ind]+h_prev[ind])/2.0)
        
    ## Calculating Electric Field and Electric Displacement
    E = -(V_grid[1:]-V_grid[:-1])/h
    D = E*perm_grid[:-1]
    
    ## Changing the units of surface densities to 1/cm2
    conc_SiO = conc_SiO*(Na/1e4)
    conc_SiOH = conc_SiOH*(Na/1e4)
    conc_SiOH2 = conc_SiOH2*(Na/1e4)

    

    return space_grid, end_s, end_o, V_grid, E, D, conc_h, conc_e, conc_H, conc_OH, conc_posI, conc_negI, conc_SiO[end_o], conc_SiOH[end_o], conc_SiOH2[end_o], rho_SC






