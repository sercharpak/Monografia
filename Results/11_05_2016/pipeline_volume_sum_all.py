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
# USAGE = "python pipeline_volume_sum_all.py OPT"
USAGE = "python pipeline_volume_sum_all.py OPT(0 or 1)"
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
if(len(sys.argv)!=2):
    print "Please use correctly"
    print USAGE
    sys.exit()

if(float(sys.argv[1])==1.0):

    array_seeds_FA = np.array([0.5,0.6,0.7,0.8,0.9])
    seed_Trace = 1.0
    #Files with the info of the search space
    #../FoF/src/grid_FA_0.5_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_0.0.dat ../FoF/src/grid_FA_0.7_Trace_0.0.dat
    #../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_1.0_Trace_0.0.dat
    #It's the first line of each file.

    #Files with the groups info
    #group_results_seeds_FA_0.5_Trace_1.0_search_FA_0.5_Trace_0.0.dat

    #Here it is done by hand

    fig = plt.figure(figsize = (8,8))
    plt.grid()
    plt.ylabel("$Volume_{max}$ / Search Space", fontsize=20)
    plt.xlabel("Growth FA Thresholds ", fontsize=20)
    plt.xlim(0.5,1.0)
    plt.ylim(0.0, 1.0)
    plt.tick_params(axis='both', which='major', labelsize=20)
    #plt.title("Group growth, search space vs search space FA threshold for \n seeds: FA: " + str(seed_FA) + " Trace: " + str(seed_Trace))
    for seed_FA in array_seeds_FA:
        print "The seed FA is", seed_FA
        results_x = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        results_y = np.zeros(6)

        i=0
        for search_FA in results_x:
            search_space_file = "../FoF/src/grid_FA_"+str(search_FA)+"_Trace_0.0.dat"
            search_space = float(readFirstLine(search_space_file))
            print "The search space for search_FA", search_FA, "is", search_space
            groups_file = "group_results_seeds_FA_"+str(seed_FA)+"_Trace_"+str(seed_Trace)+"_search_FA_"+str(search_FA)+"_Trace_0.0.dat"
            print "The first line of the group results file is:", readFirstLine(groups_file)
            group_data = np.loadtxt(groups_file)
            print "The number of groups is:", file_len(groups_file)
            if(file_len(groups_file)==1):
                #There is only one group the zeros
                group_max_vol = 0.0
                group_max_vol_number=0.0
            else:
                group_data = np.loadtxt(groups_file, skiprows=1)
                if(file_len(groups_file)>2.0):
                    volumes = group_data[:,7]
                    group_max_vol = np.max(volumes)
                    group_max_vol_ind = np.argmax(volumes)
                    group_max_vol_number = group_data[group_max_vol_ind,8]
                else:
                    volumes = group_data[7]
                    group_max_vol = volumes
                    group_max_vol_number = group_data[8]

            print "The max volume detected is: ",group_max_vol ,"for group number: ",group_max_vol_number
            results_y[i] = group_max_vol/search_space
            i+=1

        plt.plot(results_x,results_y, label="Seed FA "+str(seed_FA))

    plt.legend(loc='upper left')
    plt.savefig("./volumes_growth_FA_All_NO_ZERO.png",format = "png")
    plt.close(fig)

elif(float(sys.argv[1])==0.0):
    array_seeds_FA = np.array([0.5,0.6,0.7,0.8,0.9])
    seed_Trace = 1.0
    #Files with the info of the search space
    #../FoF/src/grid_FA_0.5_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_0.0.dat ../FoF/src/grid_FA_0.7_Trace_0.0.dat
    #../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_1.0_Trace_0.0.dat
    #It's the first line of each file.

    #Files with the groups info
    #group_results_seeds_FA_0.5_Trace_1.0_search_FA_0.5_Trace_0.0.dat

    #Here it is done by hand

    fig = plt.figure(figsize = (8,8))
    plt.grid()
    plt.ylabel("$Volume_{max}$ / Search Space", fontsize=20)
    plt.xlabel("Growth FA Thresholds ", fontsize=20)
    plt.xlim(0.5,1.0)
    plt.ylim(0.0, 1.0)
    plt.tick_params(axis='both', which='major', labelsize=20)
    #plt.title("Group growth, search space vs search space FA threshold for \n seeds: FA: " + str(seed_FA) + " Trace: " + str(seed_Trace))
    for seed_FA in array_seeds_FA:
        print "The seed FA is", seed_FA
        results_x = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        results_y = np.zeros(6)

        i=0
        for search_FA in results_x:
            search_space_file = "../FoF/src/grid_FA_"+str(search_FA)+"_Trace_0.0.dat"
            search_space = float(readFirstLine(search_space_file))
            print "The search space for search_FA", search_FA, "is", search_space
            groups_file = "group_results_seeds_FA_"+str(seed_FA)+"_Trace_"+str(seed_Trace)+"_search_FA_"+str(search_FA)+"_Trace_0.0.dat"
            print "The first line of the group results file is:", readFirstLine(groups_file)
            group_data = np.loadtxt(groups_file)
            print "The number of groups in the group results file is:", file_len(groups_file)
            if(file_len(groups_file)==1):
                #There is only one group. The 0.
                volumes = group_data[7]
                group_max_vol = volumes
                group_max_vol_number = group_data[8]
            else:
                volumes = group_data[:,7]
                group_max_vol = np.max(volumes)
                group_max_vol_ind = np.argmax(volumes)
                group_max_vol_number = group_data[group_max_vol_ind,8]

            print "The max volume detected is: ",group_max_vol ,"for group number: ",group_max_vol_number
            results_y[i] = group_max_vol/search_space
            i+=1

        plt.plot(results_x,results_y, label="Seed FA "+str(seed_FA))

    plt.legend(loc='upper left')
    plt.savefig("./volumes_growth_FA_All_WITH_ZERO.png",format = "png")
    plt.close(fig)
else:
    #Use correctly
    print "Please use correctly"
    print USAGE
    sys.exit()
