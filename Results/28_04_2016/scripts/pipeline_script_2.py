#!/usr/bin/python
#
#Written by Sergio Hernandez Charpak
#

# Usage
# python pipeline_script_2.py file_search file_seeds
#

USAGE = "python pipeline_script_2.py file_search file_seeds"

#Imports
#General Imports 
#Reading the CIC
from ast import literal_eval
from struct import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#Arguments
import sys
#Running bash
import subprocess

#Function to extract the thresholds from a filename
def get_thresholds(filename):
    file_Dat = filename.strip('.dat')
    file_array = file_Dat.split('_')
    thresh_FA = float(file_array[2])
    thresh_Trace = float(file_array[4])
    return thresh_FA, thresh_Trace

#---------------------------------------------------------

if(len(sys.argv)!=3):
    print "Please use correctly"
    print USAGE
    sys.exit() 

#Typical file_name: grid_FA_0.6_Trace_1.0.dat

#Extracts the thresholds from the filename

file_search = sys.argv[1]
file_seeds = sys.argv[2]
search_thresh_FA, search_thresh_Trace = get_thresholds(file_search)
seeds_thresh_FA, seeds_thresh_Trace = get_thresholds(file_seeds)

print "Search file thresholds (FA, Trace): ", search_thresh_FA, search_thresh_Trace
print "Seeds file thresholds (FA, Trace): ", seeds_thresh_FA, seeds_thresh_Trace

#Excecute the FoF script
# -e 0.4 means that the linking length is 0.5
# -m 10 means that it keeps groups of at least 10 particles
comando = '../FoF/src/fof -e 1.1 -m 10 < '+file_search 
process =subprocess.Popen(comando,stdout=subprocess.PIPE, stderr=None, shell=True)
resultsString=process.communicate()
print resultsString

#TODO
#loads the data from the results in the halo finder
fof_groups = loadtxt('fof.grp', skiprows=1)

#laod the mock data we wrote before
particles = loadtxt(file_search, skiprows=6)

print fof_groups.size 
print fof_groups #each entry correesponds to the ID of the fof group
id_groups = list(set(fof_groups))
print size(id_groups), 'groups'#this is the number of groups found in the FOF

index = where(fof_groups==0)
index = index[0]
x_group = particles[index,0]
y_group = particles[index,1]
z_group = particles[index,2]

#plt.figure()
#plt.scatter(x_group, z_group)
#plt.show()
