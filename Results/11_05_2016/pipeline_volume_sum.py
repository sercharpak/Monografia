#!/usr/bin/python

#Plotting the detected groups
#By Sergio Daniel Hernandez Charpak

#-------------------------------------------
#Imports
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import sys
#-------------------------------------------
# USAGE = "python pipeline_volume_sum.py 0.6 1.0 "
USAGE = "python pipeline_volume_sum.py 0.6 1.0"
#-------------------------------------------
#Functions
def readFirstLine(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
        return first_line
def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1
#-------------------------------------------
if(len(sys.argv)!=3):
    print "Please use correctly"
    print USAGE
    sys.exit()

seed_FA = float(sys.argv[1])
seed_Trace = float(sys.argv[2])

#Files with the info of the search space
#../FoF/src/grid_FA_0.5_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_0.0.dat ../FoF/src/grid_FA_0.7_Trace_0.0.dat
#../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_1.0_Trace_0.0.dat
#It's the first line of each file.

#Files with the groups info
#group_results_seeds_FA_0.5_Trace_1.0_search_FA_0.5_Trace_0.0.dat

#Here it is done by hand

results_x = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
results_y = np.zeros(6)

i=0
for search_FA in results_x:
    search_space_file = "../FoF/src/grid_FA_"+str(search_FA)+"_Trace_0.0.dat"
    search_space = float(readFirstLine(search_space_file))
    print "The search space for search_FA", search_FA, "is", search_space
    groups_file = "group_results_seeds_FA_"+str(seed_FA)+"_Trace_"+str(seed_Trace)+"_search_FA_"+str(search_FA)+"_Trace_0.0.dat"
    group_data = np.loadtxt(groups_file)
    if(search_FA==1.0):
        if(file_len(groups_file)>1.0):
            volumes = group_data[:,7]
        else:
            volumes = group_data[7]
    else:
        volumes = group_data[:,7]
    group_max_vol = np.max(volumes)
    group_max_vol_ind = np.argmax(volumes)
    print "The max volume detected is: ",group_max_vol ,"for group number: ",group_max_vol_ind
    results_y[i] = group_max_vol/search_space
    i+=1

fig = plt.figure(figsize = (8,8))
plt.plot(results_x,results_y)
plt.xlabel("search FA ")
plt.ylabel("$Volume_{max}$ / Search Space")
plt.grid()
plt.title("Group growth, search space vs search space FA threshold for \n seeds: FA: " + str(seed_FA) + " Trace: " + str(seed_Trace))
plt.savefig("./volumes_growth_FA"+str(seed_FA)+"_Trace_"+str(seed_Trace)+".png",format = "png")
plt.close(fig)
