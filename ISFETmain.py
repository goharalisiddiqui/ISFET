import ISFETmodel as isfm 
import numpy as np





L_s = 5e-7
L_o = 5e-9
L_f = 5e-7
min_spacing = 1e-11


rel_perm_s = 11.9
rel_perm_o = 3.9
rel_perm_f = 80.1

p_doping_density = 1e17                                 # Acceptor doping density of semiconductor [1/cm3]
n_doping_density = 0.0                                  # Donor doping density of semiconductor [1/cm3]


L_list = [L_s,L_o,L_f]                                  # List of lengths of the three regions
T = 300.0                                               # Absolute Temperature
rel_perm_list = [rel_perm_s,rel_perm_o,rel_perm_f]      # List of relative permittivities for the three regions                                   
salt_conc_list = [1e-3]                                 # Ionic concentrations in Molar
pH_list = 0.1*np.array(range(10,110,2))                 # pH list                                                
pK1 = 6.7                                               # pK value for SiO protonation reaction
pK2 = -1.9                                              # pK value for SiOH protonation reaction
mode = ['N', 0.0]                                       # Calculation mode options: 'N' constant bias (in V), 'B' band bending (in V) or 'C' Sheet charge (in C/cm2)
steric = True
surf_states_list = [1e14]                               # Surface sited density list (1/cm2)                     
shelf_loc = "./Results_ISFET/"                          # Location to save data 
title = "ISFET_01"                               # Name for the calculation run
force_new_calc = True                                  # Force a new calculation even if the data file already exist
plots = [True,True,True,True,True,True,True,True,True] # Plots for #elec #SSforSS #SPforSS #SpHforSS #SRforSS #SRforIC #SSforIC #SCforSS #FBforSS



isfm.isfetmodel_1D(L_list, min_spacing,rel_perm_list, salt_conc_list, pH_list, surf_states_list, p_doping_density, n_doping_density, T, pK1, pK2, mode, steric, shelf_loc, title, force_new_calc, plots)
