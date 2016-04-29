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

#Function to determine is valid or not. If it contains a seed or not
#Each time it founds a seed in a group it deletes it from the seed array to make shorter the search.
# @param int id_group the id of the group to valid
def group_valid(id_group,seeds):
    index = where(fof_groups==id_group)
    index = index[0]
    x_group = particles[index,0]
    y_group = particles[index,1]
    z_group = particles[index,2]
    index_seeds_in_group = []
    for i in range(shape(seeds)[0]):
        seed_i_x = seeds[i,0]
        seed_i_y = seeds[i,1]
        seed_i_z = seeds[i,2]
        index_seed = np.where((x_group == seed_i_x) & (y_group == seed_i_y) & (z_group == seed_i_z) )
        index_seed = index_seed[0]
        if(size(index_seed) !=0):
            #print "The seed " +str(i)+ " is in the group. Seed: ", seeds[i,:]
            index_seeds_in_group.append(i)
            #print x_group[index_seed[0]], y_group[index_seed[0]], z_group[index_seed[0]]
        #else:
            #print "The seed " +str(i)+ " is not int the group. Seed: ", seeds[i,:]
    #Now it deletes the seeds from the seed array
    if((size(index_seeds_in_group) > 0)):
        return True, index_seeds_in_group
    else:
        return False, index_seeds_in_group

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
#The search universe
particles = loadtxt(file_search, skiprows=6)
#The seeds to discriminate the groups
seeds = loadtxt(file_seeds, skiprows=6)

print fof_groups.size
print fof_groups #each entry correesponds to the ID of the fof group
id_groups = list(set(fof_groups))
print size(id_groups), 'groups'#this is the number of groups found in the FOF

#output_name = './fof_seeds_FA_'+str(seeds_thresh_FA)+'_Trace_'+str(seeds_thresh_Trace)+'_search_FA_'str(search_thresh_FA)+'_Trace_'+str(search_thresh_Trace)+'.dat'
#print "Writing in file " + output_name
#fileout = open(output_name, 'w')
#--------------------------------------
#----------------------------------------
valid_id_groups = []

print shape(seeds), 'shape seeds'
for id_group in id_groups:
    print "There are:", shape(seeds)[0], "seeds left to analyze"
    valid_condition, index_to_delete = group_valid(id_group,seeds)
    print "Ended the evaluation for group number", id_group
    if(valid_condition):
        #It is a valid group
        valid_id_groups.append(id_group)
        print "number of seeds to delete", len(index_to_delete)
        seeds = np.delete(seeds, index_to_delete, axis=0)
    #else:
        #It is not a valid group
    if(shape(seeds)[0]==0):
        #No more seeds to evaluate. Can break
        break
print "There are ", size(valid_id_groups), "valid groups"
print "Oposite to the ", size(id_groups), "groups at first"



#plt.figure()
#plt.scatter(x_group, z_group)
#plt.show()
