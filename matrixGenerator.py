#! /usr/bin/python3
import numpy as np
from scipy.sparse import diags




def get_capMat_1D_spherical(space_grid, perm_grid, BC_left, BC_right):
    
    space_grid = np.array(space_grid)
    perm_grid = np.array(perm_grid)

    h = np.append(space_grid[1:] - space_grid[:-1],[0.0])
    h_prev = np.append([0.0],h[:-1])
    h_avg = (h+h_prev)/2.0

    p_par = (np.power(space_grid+(h/2),3) - np.power(space_grid-(h_prev/2),3))/3.0


    ## Main Diagonal
    Di0_1 = np.power(space_grid-(h_prev/2),2)
    Di0_2 = np.power(space_grid+(h/2),2)
    Di0_3 = Di0_1 + Di0_2
    Di0_4 = (Di0_3*perm_grid)/h_avg
    DiM   = np.divide(Di0_4,p_par)

                                                                                                    

    ## Right diagonal
    DiR = -np.divide(((perm_grid/h_avg)*Di0_2),p_par)[:-1]



    ## right diagonal BC
    if BC_left == 'D':
        DiR[0] = 0.0                                                            #Dirichlet for the left boundary
        DiM[0] = 1.0                                                            



    ## Left diagonal
    DiL = -np.divide(((perm_grid/h_avg)*Di0_1),p_par)[1:] 



    ## left diagonal BC
    if BC_right == 'D':
        DiL[-1] = 0.0                                                          #Dirichlet condition for right boundary
        DiM[-1] = 1.0

    cap_mat = diags([DiL,DiM,DiR], [-1,0,1], shape=(space_grid.size, space_grid.size))

    return cap_mat


def get_capMat_1D_cartesian(space_grid, perm_grid, BC_left, BC_right):
    
    ## Spatial coordinates of the grid points
    space_grid = np.array(space_grid)

    
    ## Distance of each grid point with its right neighbour, zero for last point to make it same size as space grid
    h = np.append(space_grid[1:] - space_grid[:-1],[0.0])
    ## Distance of each grid point with its left neighbour, zero for first point to make it same size as space_grid
    h_prev = np.append([0.0],h[:-1])
    
    ## Finite volume of each grid point
    h_avg = (h+h_prev)/2.0
    
    
    ## Permittivity of space between each grid point and its right neighbour
    perm_grid = np.array(perm_grid)
    ## Permittivity of space between each grid point and its left neighbour
    perm_grid_prev = np.append([0.0],perm_grid[:-1])


    ## Main Diagonal
    DiM = (1.0/h_avg[1:-1])* ((perm_grid_prev[1:-1]/h_prev[1:-1]) + (perm_grid[1:-1]/h[1:-1]))
    # First and Last grid points with only one neighbour
    DiM = np.append([(1.0/h_avg[0])*(perm_grid[0]/h[0])] , np.append(DiM , [(1.0/h_avg[-1])*(perm_grid_prev[-1]/h_prev[-1])]))

                                                                                                    

    ## Right diagonal
    DiR = -(1.0/h_avg[:-1]) * (perm_grid[:-1]/h[:-1])



    ## right diagonal BC (left BC)
    if BC_left == 'D':
        DiR[0] = 0.0                                                            #Dirichlet for the left boundary
        DiM[0] = 1.0                                                            



    ## Left diagonal
    DiL = -(1.0/h_avg[1:]) * (perm_grid_prev[1:]/h_prev[1:]) 



    ## left diagonal BC (right BC)
    if BC_right == 'D':
        DiL[-1] = 0.0                                                          #Dirichlet condition for right boundary
        DiM[-1] = 1.0

    ## Sparse capacitance matrix
    cap_mat = diags([DiL,DiM,DiR], [-1,0,1])

    return cap_mat
