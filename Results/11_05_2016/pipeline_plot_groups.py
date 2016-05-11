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
# USAGE = "python pipepline_plt_groups.py file_valid_groups output_folder OPT "
#File example: group_valid_positions_FA_Trace_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat
USAGE = "python pipepline_plt_groups.py file_valid_groups output_folder OPT(1 for plotting all 0 for not) "
SCALE_FACTOR = 976.5625
GRID = 256
#-------------------------------------------
#Functions
#Function to extract the thresholds from a group filename
def get_thresholds(filename):
    file_Dat = filename.strip(".dat")
    file_array = file_Dat.split("_")
    thresh_search_FA = float(file_array[12])
    thresh_search_Trace = float(file_array[14])
    thresh_seeds_FA = float(file_array[7])
    thresh_seeds_Trace = float(file_array[9])
    return thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace
#-------------------------------------------
if(len(sys.argv)!=4):
    print "Please use correctly"
    print USAGE
    sys.exit()

inputfile = sys.argv[1]
output_folder = sys.argv[2]
print "Processing the file: " + inputfile
thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(inputfile)
group_data = np.loadtxt(inputfile)
n_particles = len(group_data)
print "There are", n_particles,"particles for the file", inputfile

print "The grid size is: ", GRID
print "The scale factor is: ", SCALE_FACTOR
print "The real boxsize is [kpc/h]:", GRID*SCALE_FACTOR, " or in [Mpc/h]: ", (GRID*SCALE_FACTOR)/1000.0

#------------------------------------
#First we plot all the particles
from mpl_toolkits.mplot3d import Axes3D
#Particles coordinates
x_particles_all = group_data[:,0]
y_particles_all = group_data[:,1]
z_particles_all = group_data[:,2]
#We transform the coordinates from GRID to the real units in Mpc/h
x_particles_all = (x_particles_all*SCALE_FACTOR)/1000.0
y_particles_all = (y_particles_all*SCALE_FACTOR)/1000.0
z_particles_all = (z_particles_all*SCALE_FACTOR)/1000.0
group_particle = group_data[:,5]
n_groups = len(group_particle)
#TODO This should be improved
n_x, n_y, n_z = (GRID*SCALE_FACTOR)/1000.0,(GRID*SCALE_FACTOR)/1000.0,(GRID*SCALE_FACTOR)/1000.0

#This is the plot with the noise. The 0 group is the one with the particles which do not belong to a group
"""
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_particles_all,y_particles_all,z_particles_all,s=2, c=group_particle)
ax.set_xlim3d(0, n_x)
ax.set_ylim3d(0, n_y)
ax.set_zlim3d(0, n_z)
ax.set_xlabel('x [Mpc/h]', size=30)
ax.set_ylabel('y [Mpc/h]', size=30)
ax.set_zlabel('z [Mpc/h]', size=30)
ax.set_title('Regions \n Seeds (FA:'+str(thresh_seeds_FA)+
' Divergence:'+str(thresh_seeds_Trace)+') \n Growth (FA:'+str(thresh_search_FA)+
' Divergence:'+str(thresh_search_Trace)+')', size=30)
ax.tick_params(axis='both', which='major', labelsize=25)
plt.savefig(output_folder+"regions_all_seeds_FA_" + str(thresh_seeds_FA) + "_Trace_" + str(thresh_seeds_Trace) + "_search_FA_" + str(thresh_search_FA) + "_Trace_" + str(thresh_search_Trace) + ".png",format = 'png')
plt.close(fig)
"""
#The group 0.0 seems like noise. We will proceed to make the 3d plot of the regions excluding this one.
index = where(group_particle==0.0)
index = index[0]
x_particles_all = np.delete(x_particles_all, index, axis=0)
y_particles_all = np.delete(y_particles_all, index, axis=0)
z_particles_all = np.delete(z_particles_all, index, axis=0)
group_particle = np.delete(group_particle, index, axis=0)
#The plot
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_particles_all,y_particles_all,z_particles_all,s=2, c=group_particle)
ax.set_xlim3d(0, n_x)
ax.set_ylim3d(0, n_y)
ax.set_zlim3d(0, n_z)
ax.set_xlabel('x [Mpc/h]', size=30)
ax.set_ylabel('y [Mpc/h]', size=30)
ax.set_zlabel('z [Mpc/h]', size=30)
ax.set_title('Regions \n Seeds (FA:'+str(thresh_seeds_FA)+
' Divergence:'+str(thresh_seeds_Trace)+') \n Growth (FA:'+str(thresh_search_FA)+
' Divergence:'+str(thresh_search_Trace)+')', size=30)
ax.tick_params(axis='both', which='major', labelsize=25)
plt.savefig(output_folder+"regions_nonoise_seeds_FA_" + str(thresh_seeds_FA) + "_Trace_" + str(thresh_seeds_Trace) + "_search_FA_" + str(thresh_search_FA) + "_Trace_" + str(thresh_search_Trace) + ".png",format = 'png')
plt.close(fig)
#-------------------------------
#If the option 1 was selected we plot each of the regions
#Now let us plot each region.
print "The option selected is", sys.argv[3]
if(int(sys.argv[3]) == 1):
    id_groups = list(set(group_particle))
    id_groups = np.array(id_groups)
    for id_group in id_groups:
        print "Is processing group number", id_group
        index = where(group_particle==id_group)
        index = index[0]
        x_particles = x_particles_all[index]
        y_particles = y_particles_all[index]
        z_particles = z_particles_all[index]
        x_max = np.max(x_particles)
        x_min = np.min(x_particles)
        y_max = np.max(y_particles)
        y_min = np.min(y_particles)
        z_max = np.max(z_particles)
        z_min = np.min(z_particles)
        #Only the region in its region
        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_particles,y_particles,z_particles,s=15)
        ax.set_xlim3d(x_min, x_max)
        ax.set_ylim3d(y_min, y_max)
        ax.set_zlim3d(z_min, z_max)
        ax.set_xlabel('x [Mpc/h]', size=30)
        ax.set_ylabel('y [Mpc/h]', size=30)
        ax.set_zlabel('z [Mpc/h]', size=30)
        ax.set_title('Region '+str(id_group)+' \n Seeds (FA:'+str(thresh_seeds_FA)+
        ' Divergence:'+str(thresh_seeds_Trace)+') \n Growth (FA:'+str(thresh_search_FA)+
        ' Divergence:'+str(thresh_search_Trace)+')', size=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.savefig(output_folder+"region_"+str(id_group)+"_small_seeds_FA_" + str(thresh_seeds_FA) + "_Trace_" + str(thresh_seeds_Trace) + "_search_FA_" + str(thresh_search_FA) + "_Trace_" + str(thresh_search_Trace) + ".png",format = 'png')
        plt.close(fig)
        #At scale with the simulation
        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_particles,y_particles,z_particles,s=2)
        ax.set_xlim3d(0, n_x)
        ax.set_ylim3d(0, n_y)
        ax.set_zlim3d(0, n_z)
        ax.set_xlabel('x [Mpc/h]', size=30)
        ax.set_ylabel('y [Mpc/h]', size=30)
        ax.set_zlabel('z [Mpc/h]', size=30)
        ax.set_title('Region '+str(id_group)+' at scale \n Seeds (FA:'+str(thresh_seeds_FA)+
        ' Divergence:'+str(thresh_seeds_Trace)+') \n Growth (FA:'+str(thresh_search_FA)+
        ' Divergence:'+str(thresh_search_Trace)+')', size=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        plt.savefig(output_folder+"region_"+str(id_group)+"_scale_seeds_FA_" + str(thresh_seeds_FA) + "_Trace_" + str(thresh_seeds_Trace) + "_search_FA_" + str(thresh_search_FA) + "_Trace_" + str(thresh_search_Trace) + ".png",format = 'png')
        plt.close(fig)
