#!/usr/bin/python

#Plotting the general properties of the detected groups
#By Sergio Daniel Hernandez Charpak

#-------------------------------------------
#Imports
import numpy as np
import matplotlib.pyplot as plt
import sys
#-------------------------------------------
# Constants
# USAGE = "python pipeline_histograms.py file_groups"
USAGE = "python pipeline_histograms.py file_groups"
SCALE_FACTOR = 976.5625
GRID = 256
VOLUME = 256.0 ** 3.0


#Function to extract the thresholds from a group filename
def get_thresholds(filename):
    file_Dat = filename.strip(".dat")
    file_array = file_Dat.split("_")
    thresh_search_FA = float(file_array[9])
    thresh_search_Trace = float(file_array[11])
    thresh_seeds_FA = float(file_array[4])
    thresh_seeds_Trace = float(file_array[6])
    return thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace

def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

MPC_SCALE_FACTOR =  (SCALE_FACTOR)/1000.0
#Seed 6, 7, 8, 9 Growth 9
array_seeds_FA = np.array([0.6,0.7,0.8,0.9])
#group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat


#with Laniakea
#Radius Laniakea: 80 Mpc/h
#Volume Laniakea: 2.14E16(Mpc/h)^3

radius_laniakea_Mpc = 80.0
volume_laniakea_Mpc = 2.14 * (10.0**6.0)
log_10_vol_laniakea_Mpc = np.log10(volume_laniakea_Mpc)

volume_sim_mpc = 250.0**3.0
log_10_volume_sim_mpc = np.log10(volume_sim_mpc)

fig = plt.figure(figsize = (14,14))
plt.axvline(log_10_vol_laniakea_Mpc, linewidth=4, color='m', label='Laniakea')
plt.axvline(log_10_volume_sim_mpc, linewidth=4, color='chocolate', label='Volume Simulation')


for thresh_seeds_FA in array_seeds_FA:

    inputfile = "group_results_seeds_FA_"+str(thresh_seeds_FA)+"_Trace_1.0_search_FA_0.9_Trace_0.0.dat"
    print "Processing the file: " + inputfile
    thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(inputfile)
    group_data = np.loadtxt(inputfile)
    n_groups = file_len(inputfile)
    print "There are", n_groups,"groups"


#-------------------------------------------
#Volume distribution
    if(n_groups==1):
        volumes = group_data[7]
    else:
        volumes = group_data[:,7]

    volumes = volumes * (MPC_SCALE_FACTOR**3.0)
    #print "----volumes-----"
    #print (np.log10(volumes))
    #print "----laniakea------"
    #print log_10_vol_laniakea_Mpc
    #print "----------------"

    bins=25
    #plt.hist(np.log10(volumes), bins=bins, histtype='step',linewidth=3,label = "FA = "+str(thresh_seeds_FA))

    hist_vol, bins = np.histogram(np.log10(volumes),bins=bins)

    print bins
    n_bins = np.zeros(len(hist_vol))
    for i in range(len(hist_vol)):
        n_bins[i] = (bins[i]+bins[i+1])/2.0
    #hist_array = np.array(hist_vol)
    print hist_vol
    log_10_hist_vol = np.log10(hist_vol)
    print log_10_hist_vol
    where_are_NaNs = np.isinf(log_10_hist_vol)
    log_10_hist_vol[where_are_NaNs] = 0
    print log_10_hist_vol

    plt.step( n_bins,log_10_hist_vol,linewidth=3, label = "FA = "+str(thresh_seeds_FA) )


plt.xlabel("$log_{10}$ of Volume $(Mpc/h)^3$", fontsize=20)
plt.ylabel("$log_{10}$ of Number of Groups", fontsize=20)
#plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.legend(loc='upper right', fontsize='20')
plt.title("Distribution of volume $(Mpc/h)^3$", fontsize=20)
plt.savefig("./volumes_distr_Mpc_laniakea_all_plot.png",format = "png")
#plt.savefig("./volumes_distr_Mpc_laniakea_all_histstep.png",format = "png")
plt.show()
plt.close(fig)
