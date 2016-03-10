from pygadgetreader import *
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

snap_file_sim='/lustre/home/ciencias/fisica/pregrado/sd.hernandez204/Outputs/28_02_2016/snapshot_005'

requested='dmcount'
n_dm=readheader(snap_file_sim,requested)
print "Number of Dark Matter particles"
print n_dm

array_pos=readsnap(snap_file_sim,'pos','dm')

array_vel=readsnap(snap_file_sim,'vel','dm')




array_magnitud_vel=np.zeros(n_dm)
for i in range(0,n_dm):
    array_magnitud_vel[i]=np.linalg.norm(array_vel[i,:])
print "The speed array has been calculated"

longitud = (int) (array_pos.shape[0])
print "The lenght of the position array should be the same number of particles"
print str(longitud)


array_copy=array_pos.copy()
array_results=np.zeros((longitud,4))
array_results[:,:-1] = array_copy
print "The results array has been created."


vel_max = np.max(array_magnitud_vel)
vel_min = np.min(array_magnitud_vel)
vel_std = np.std(array_magnitud_vel)
vel_mean = np.mean(array_magnitud_vel)

print "Vel_max, Vel_min, std_vel, vel_mean"
print vel_max, vel_min, vel_std, vel_mean
thresh_v_hig = vel_mean + 2.0 * vel_std
thresh_v_low = vel_mean -  vel_std
print "Thresh_high, Thresh_low"
print thresh_v_hig, thresh_v_low

#Gets the indexes of the seeds
indexes = np.where(array_magnitud_vel>(thresh_v_hig))[0]
print 'There are: '+ str(indexes.shape[0]) +' seeds'


def region_accretion(l):
    #Resize the window if necesary
    size_window_up = n_window
    if(((l+size_window_up) > longitud)):
        size_window_up = longitud - l
    size_window_down = n_window
    if(((l-size_window_down) < 0)):
        size_window_down = 0 + l
    distancias_l_up = np.zeros(size_window_up)
    distancias_l_down = np.zeros(size_window_down)
    distancias_l_up[:] = 1000000
    distancias_l_down[:] = 1000000
    if(size_window_up == size_window_down):
        for k in range(0,size_window_up):
            if(((l+k) < longitud)):
                if((array_results[l+k,3] == 0)):
                    distancias_l_up[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l+k,0:2])**(2.0) ))
            if(((l-k) > 0)):
                if((array_results[l-k,3] == 0)):
                    distancias_l_down[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l-k,0:2])**(2.0) )) 
    else:
        for k in range(0,size_window_up):
            if(((l+k) < longitud)):
                if((array_results[l+k,3] == 0)):
                    distancias_l_up[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l+k,0:2])**(2.0) ))
        for k in range(0,size_window_up):
            if(((l-k) > 0)):
                if((array_results[l-k,3] == 0)):
                    distancias_l_down[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l-k,0:2])**(2.0) ))
    distancias_l_up[0]=1000000
    distancias_l_down[0]=1000000
    closest_l_up = np.min(distancias_l_up)
    index_closest_up = np.argmin(distancias_l_up)
    closest_l_down = np.min(distancias_l_down)
    index_closest_down = np.argmin(distancias_l_down)
    #Now si if these validate the condition 
    #It has to see if it got any new candidates (not marked particles)
    if ((closest_l_up != 1000000)):
        vel_l = array_magnitud_vel[l]
        vel_closest_up = array_magnitud_vel[l + index_closest_up]
        if( (vel_l >= vel_closest_up >= thresh_v_low) & (array_results[l + index_closest_up,3] == 0) ):
            #it is marked
            array_results[l + index_closest_up,3] = 1
            #Recursion
            l_up = l + index_closest_up
            region_accretion(l_up)
    if ((closest_l_down != 1000000)):
        vel_l = array_magnitud_vel[l]
        vel_closest_down = array_magnitud_vel[l + index_closest_down]
        if( (vel_l >= vel_closest_down >= thresh_v_low) & (array_results[l - index_closest_down,3] == 0) ):
            #it is marked
            array_results[l - index_closest_down,3] = 1
            #Recursion
            l_down = l - index_closest_down
            region_accretion(l_down)


#Define the default size of the search window
n_window=10000
#Define the minimum distance
#TODO CHANGE THIS
#minimum_distance=1500

#Now mark the seeds
for i in indexes:
    array_results[i,3] = 1
#Now do the recursion
for i in indexes:
    #Recursion
    region_accretion(i)

print 'The regions have: '+ str(len(array_results[array_results[:,3] == 1,:])) +' galaxies.'
print 'Starting with points with V > '+ str(thresh_v_hig)
print 'Ending at points with V < ' + str(thresh_v_low)

#Now it gets the points of the regions
array_pos_region=array_results[(array_results[:,3] == 1),:]

#It gets the statistics out of the array
results_mean = np.mean(array_pos_region)
results_std = np.std(array_pos_region)
print "mean, std"
print results_mean, results_std
#Now the statistics out of the velocities
array_vel_region=array_magnitud_vel[(array_results[:,3] == 1)]
results_vel_mean = np.mean(array_vel_region)
results_vel_std = np.std(array_vel_region)
print "vel_mean, vel_std"
print results_vel_mean, results_vel_std

#Resize of the box
l_x_max = np.max(array_pos_region[:,0])
l_y_max = np.max(array_pos_region[:,1])
l_z_max = np.max(array_pos_region[:,2])
l_x_min = np.min(array_pos_region[:,0])
l_y_min = np.min(array_pos_region[:,1])
l_z_min = np.min(array_pos_region[:,2])

#Cube which has the region
line_x = [l_x_min,l_x_max, l_x_max, l_x_min, l_x_min, l_x_min, l_x_min, l_x_min, l_x_min, l_x_max, l_x_max, l_x_max, l_x_max, l_x_max, l_x_max, l_x_min]
line_y = [l_y_min,l_y_min, l_y_max, l_y_max, l_y_min, l_y_min, l_y_max, l_y_max, l_y_max, l_y_max, l_y_max, l_y_max, l_y_min, l_y_min, l_y_min, l_y_min]
line_z = [l_z_min,l_z_min, l_z_min, l_z_min, l_z_min, l_z_max, l_z_max, l_z_min, l_z_max, l_z_max, l_z_min, l_z_max, l_z_max, l_z_min, l_z_max, l_z_max]

fig = plt.figure(figsize=(25, 25))
ax = fig.add_subplot(111, projection='3d')
ax.plot(line_x, line_y, line_z, color = 'r')
ax.scatter(array_pos_region[:,0],array_pos_region[:,1],array_pos_region[:,2],s=20,c='b')
ax.set_xlim3d(l_x_min, l_x_max)
ax.set_ylim3d(l_y_min,l_y_max)
ax.set_zlim3d(l_z_min,l_z_max)
ax.set_xlabel('X kpc/h', size=30)
ax.set_ylabel('Y kpc/h', size=30)
ax.set_zlabel('Z kpc/h', size=30)
ax.set_title('Regions obtained', size=35)
ax.tick_params(axis='both', which='major', labelsize=30)
plt.savefig("./regions_small.png",format = 'png')
plt.close(fig)

fig = plt.figure(figsize=(25, 25))
ax = fig.add_subplot(111, projection='3d')
ax.plot(line_x, line_y, line_z, color = 'r')
ax.scatter(array_pos_region[:,0],array_pos_region[:,1],array_pos_region[:,2],s=20)
ax.set_xlim3d(0, 500000)
ax.set_ylim3d(0, 500000)
ax.set_zlim3d(0, 500000)
ax.set_xlabel('X kpc/h', size=30)
ax.set_ylabel('Y kpc/h', size=30)
ax.set_zlabel('Z kpc/h', size=30)
ax.set_title('Regions obtained at scale with the simulation', size=35)
ax.tick_params(axis='both', which='major', labelsize=30)
#plt.colorbar()
#plt.show()
plt.savefig("./regions_big.png",format = 'png')
plt.close(fig)

#Array of points not in the region in the cube.
array_pos_plot_region = array_results[((array_results[:,0]>l_x_min) & (array_results[:,0]<l_x_max) & (array_results[:,1]>l_y_min) & (array_results[:,1]<l_y_max) & (array_results[:,2]<l_z_max) & (array_results[:,2]>l_z_min) & (array_results[:,3] == 0) ),:]

fig = plt.figure(figsize=(25, 25))
ax = fig.add_subplot(111, projection='3d')
ax.plot(line_x, line_y, line_z, color = 'r')
ax.scatter(array_pos_plot_region[:,0],array_pos_plot_region[:,1],array_pos_plot_region[:,2],s=2,c='r')
ax.scatter(array_pos_region[:,0],array_pos_region[:,1],array_pos_region[:,2],s=20, c='b')
ax.set_xlim3d(l_x_min, l_x_max)
ax.set_ylim3d(l_y_min,l_y_max)
ax.set_zlim3d(l_z_min,l_z_max)
ax.set_xlabel('X kpc/h', size=30)
ax.set_ylabel('Y kpc/h', size=30)
ax.set_zlabel('Z kpc/h', size=30)
ax.set_title('Regions obtained at scale, \n with the other points of the simulation in the region', size=35)
ax.tick_params(axis='both', which='major', labelsize=30)
plt.savefig("./regions_all_small.png",format = 'png')
plt.close(fig)

del array_pos_plot_region

#Array of points not in the region in the simulation
array_pos_plot_region = array_results[(array_results[:,3] == 0 ),:]

fig = plt.figure(figsize=(30, 30))
ax = fig.add_subplot(111, projection='3d')
ax.plot(line_x, line_y, line_z, color = 'r')
ax.scatter(array_pos_plot_region[:,0],array_pos_plot_region[:,1],array_pos_plot_region[:,2],s=1,c='r')
ax.scatter(array_pos_region[:,0],array_pos_region[:,1],array_pos_region[:,2],s=20, c='b')
ax.set_xlim3d(0, 500000)
ax.set_ylim3d(0, 500000)
ax.set_zlim3d(0, 500000)
ax.set_xlabel('X kpc/h', size=30)
ax.set_ylabel('Y kpc/h', size=30)
ax.set_zlabel('Z kpc/h', size=30)
ax.set_title('Regions obtained with the rest of the simulation points', size=35)
ax.tick_params(axis='both', which='major', labelsize=30)
plt.savefig("./regions_all_big.png",format = 'png')
plt.close(fig)
