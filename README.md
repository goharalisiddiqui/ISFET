# ISFET Model
A 1D solver for Modified-Poisson-Boltzmann Model for Ion-Sensitive Field Effect Transistor (ISFET)

# DEPENDENCIES:

These python modules must be available in your system to use this library

1. python3-numpy

2. python3-matplotlib

3. python3-scipy

4. python3-progressbar

Add the repository to your $PYTHONPATH with

$ export PYTHONPATH=$PYTHONPATH:/path/to/ISFET/

and also add this line to your ~/.bashrc or whichever shell you are using.



# USER GUIDE:
To use the model, the only file you have to edit is the ISFETmain.py file where you specify all the calculation parameters described below including the location to save the output plots.

### Inputs:

#### Required:

L_s = Length in m of the semiconductor region

L_o = Length in m of the oxide region

L_f = Length in m of the fluid region

min_spacing = Minimum grid spacing(resolution) in m used to make the spatial grid

rel_perm_s = Relative permittivity of the semiconductor region

rel_perm_o = Relative permittivity of the oxide region

rel_perm_f = Relative permittivity of the fluid region

p_doping_density = Acceptor doping density of semiconductor in 1/cm3

n_doping_density = Donor doping density of semiconductor in 1/cm3

salt_conc_list = List of Ionic concentrations in Molar

pH_list = pH list                                                

surf_states_list = Surface sites density list in 1/cm2                     


#### Optional:

T = 300.0 = Absolute Temperature

pK1 = 6.7 = pK value for SiO protonation reaction

pK2 = -1.9 = pK value for SiOH protonation reaction

mode = ['N', 0.0] = Calculation mode options: 'N' with constant bias (in V), 'B' with band bending (in V) or 'C' with Sheet charge density (in C/cm2)

steric = True = Weather to used size-modified poisson-boltmann or not

shelf_loc = "./" = Location to save data 

title = "Untitled" = Name for the calculation run

force_new_calc = False = Force a new calculation even if the data file already exist

plots = [True,True,True,True,True,True,True,True,True] = Flags for output plots relatively for 
								1) Electrostatics for each pH
 								2) Surface States density for the three surface states for difffernet surface sites density 
								3) Surface potential different in reference to the applied fluid bias for different surface sites density
								4) Surface pH vs Bulk pH for different surface sites density
								5) Sensing ratio vs bulk pH for different surface sites density
								6) Sensing ratio vs bulk pH for differnet salt concentrations 
 								7) Surface States density for the three surface states for difffernet salt concentrations 
								8) Sheet charge density vs pH for different surface sites density  
								9) Fluid bias applied vs pH for different surface sites density

    
### Outputs: 

Plots the specified data plots and saves them in the specified location along with a .out file containing all the data of the calculation run


## Warnings:

** The arguments are not yet robust because there is no checking of inputs implemented yet so a wrong or unphysical values may still produce results so check the input values. **

