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
steric = True
surf_states_list = [1e13]                               # Surface sited density list (1/cm2)                     
shelf_loc = "./Results_ISFET/"                          # Location to save data 
force_new_calc = True                                  # Force a new calculation even if the data file already exist
plots = [False,False,False,False,False,False,False,False,False] # Plots for #elec #SSforSS #SPforSS #SpHforSS #SRforSS #SRforIC #SSforIC #SCforSS #FBforSS


charge_values = np.linspace(0.0,-1e-6,20)
for ind,charge in enumerate(charge_values):
    mode = ['C',charge]
    if ind <= 9:
        title = "ISFET_IV_00"+str(ind)
    else:
        title = "ISFET_IV_0"+str(ind)
    isfm.isfetmodel_1D(L_list, min_spacing,rel_perm_list, salt_conc_list, pH_list, surf_states_list, p_doping_density, n_doping_density, T, pK1, pK2, mode, steric, shelf_loc, title, force_new_calc, plots)
    print("\n\n Done with caculation %d out of %d \n\n"%(ind+1,charge_values.size))
