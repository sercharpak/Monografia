from ast import literal_eval
from struct import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# A Region Accretion Algorithm using the FA and the VDH

#Predicate to be part of the group
def predicate(i,j,k, seed_FA):
    if(results[i,j,k]>=1):
        return False #There is no need to re-evaluate the cell
    elif( (FA[i,j,k] <= thresh_FA) & (Trace[i,j,k] >= thresh_Trace)):
    #elif((Trace[i,j,k] >= thresh_Trace)):
        #It validates the predicate
        return True
    else:
        #It fails the predicate
        return False

#The algorithm, it evaluates the special cases
def region_accretion_4_neigboors(i,j,k):
    present_FA = FA[i,j,k]
    present_region = results[i,j,k]
    #Finds if it is on the border of the grid
    case = 0
    #Cases
    if(i==0):
        case+=1
    if(j==0):
        case+=2
    if(k==0):
        case+=4
    if(i==n_grid):
        case+=8
    if(j==n_grid):
        case+=16
    if(k==n_grid):
        case+=32
    #General case
    if(case==0):
        if(predicate(i+1,j,k,present_FA)):
            results[i+1,j,k] = present_region
            region_accretion_4_neigboors(i+1,j,k)
        if(predicate(i,j+1,k,present_FA)):
            results[i,j+1,k] = present_region
            region_accretion_4_neigboors(i,j+1,k)
        if(predicate(i,j,k+1,present_FA)):
            results[i,j,k+1] = present_region
            region_accretion_4_neigboors(i,j,k+1)
        if(predicate(i,j-1,k,present_FA)):
            results[i,j-1,k] = present_region
            region_accretion_4_neigboors(i,j-1,k)
        if(predicate(i,j,k-1,present_FA)):
            results[i,j,k-1] = present_region
            region_accretion_4_neigboors(i,j,k-1)
        if(predicate(i-1,j,k,present_FA)):
            results[i-1,j,k] = present_region
            region_accretion_4_neigboors(i-1,j,k)
    #Special cases
    #... Border cases

inputfolder = 'simulation_folder/'
inputfile_1 = 'snapshot_005.eigen_1'
inputfile_2 = 'snapshot_005.eigen_2'
inputfile_3 = 'snapshot_005.eigen_3'

# Formation of the Shear Tensor eigenvalues grids, FA and VDH
grid_1 = read_CIC_scalar(inputfolder+inputfile_1)
grid_2 = read_CIC_scalar(inputfolder+inputfile_2)
grid_3 = read_CIC_scalar(inputfolder+inputfile_3)

FA = (grid_1-grid_3)**2  + (grid_2-grid_3)**2  + (grid_1-grid_2)**2 
FA = FA/(grid_1**2 + grid_2**2 + grid_3**2)
FA = np.sqrt(FA)/np.sqrt(3.0)

Trace = grid_1 + grid_2 + grid_3

#Seed selection
n_x,n_y,n_z = shape(FA)
n_grid = n_x-1 # the size of box -1, the maximum index
results = np.zeros(shape(FA)) #preps the results matrix
thresh_FA = 0.8
thresh_Trace = 0.0
thresh_FA_seeds = 0.6
thresh_Trace_seeds = 1.0

#Needs to pick the seeds. 
seed_list = []

counter = 0
for i in range (n_x):
    for j in range (n_y): 
        for k in range (n_z):
            if((FA[i,j,k]<thresh_FA_seeds ) & (Trace[i,j,k]>thresh_Trace_seeds)):
                counter +=1
                results[i,j,k] = counter
                seed_list.append([i,j,k])


#Growth thresholds selection
thresh_FA = 0.8
thresh_Trace = 0.0

#Runs the Algorithm
for (i,j,k) in seed_list:
    region_accretion_4_neigboors(i,j,k)


