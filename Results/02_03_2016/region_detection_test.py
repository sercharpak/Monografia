
#This module is from https://bitbucket.org/rthompson/pygadgetreader
from pygadgetreader import *




from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np



snap_file_sim='../../Outputs/28_02_2016/snapshot_005'



requested='dmcount'
n_dm=readheader(snap_file_sim,requested)
#n_dm=readheader(snap_file_example,requested)
print n_dm



array_pos=readsnap(snap_file_sim,'pos','dm')



longitud = (int) (array_pos.shape[0])



array_vel=readsnap(snap_file_sim,'vel','dm')



array_magnitud_vel=np.zeros(n_dm)
for i in range(0,n_dm):
    array_magnitud_vel[i]=np.linalg.norm(array_vel[i,:])



array_copy=array_pos.copy()



array_results=np.zeros((longitud,4))



array_results[:,:-1] = array_copy



vel_max = np.max(array_magnitud_vel)
vel_min = np.min(array_magnitud_vel)
vel_std = np.std(array_magnitud_vel)
print vel_max, vel_min, vel_std
thresh_v_hig = vel_max-vel_std/1.0
thresh_v_low = vel_min+vel_std/5.0
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
    distancias_i_up = np.zeros(size_window_up)
    distancias_i_down = np.zeros(size_window_down)
    distancias_i_up[:] = 100000
    distancias_i_down[:] = 100000
    if(size_window_up == size_window_down):
        for k in range(0,size_window_up):
            if(((l+k) < longitud)):
                if((array_results[l+k,3] == 0)):
                    distancias_i_up[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l+k,0:2])**(2.0) ))
            if(((l-k) > 0)):
                if((array_results[l-k,3] == 0)):
                    distancias_i_down[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l-k,0:2])**(2.0) )) 
    else:
        for k in range(0,size_window_up):
            if(((l+k) < longitud)):
                if((array_results[l+k,3] == 0)):
                    distancias_i_up[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l+k,0:2])**(2.0) ))
        for k in range(0,size_window_up):
            if(((l-k) > 0)):
                if((array_results[l-k,3] == 0)):
                    distancias_i_down[k]= np.sqrt( np.sum( (array_results[l,0:2] - array_results[l-k,0:2])**(2.0) ))
    distancias_i_up[0]=100000
    distancias_i_down[0]=100000
    closest_i_up = np.min(distancias_i_up)
    index_closest_up = np.argmin(distancias_i_up)
    closest_i_down = np.min(distancias_i_down)
    index_closest_down = np.argmin(distancias_i_down)
    #print closest_i_up, closest_i_down
    #print index_closest_up, index_closest_down
    #print l
    #Now si if these validate the condition
    vel_i = array_magnitud_vel[l]
    vel_closest_up = array_magnitud_vel[l + index_closest_up]
    vel_closest_down = array_magnitud_vel[l + index_closest_down]
    if( (vel_i > vel_closest_up > thresh_v_low) & (array_results[l + index_closest_up,3] == 0) ):
        #it is marked
        array_results[l + index_closest_up,3] = 1
        #Recurtion
        region_accretion(l + index_closest_up)
    if( (vel_i > vel_closest_down > thresh_v_low) & (array_results[l - index_closest_down,3] == 0) ):
        #it is marked
        array_results[l - index_closest_down,3] = 1
        #Recurtion
        region_accretion(l - index_closest_down)




#Define the default size of the search window
n_window=10000
#Define the minimum distance
#TODO CHANGE THIS
minimum_distance=1500


#Now mark the points
for i in indexes:
    array_results[i,3]=1
    #Recursion
    region_accretion(i)




print 'The regions have: '+ str(len(array_results[array_results[:,3] == 1,:])) +' galaxies.'
print 'Starting with points with V > '+ str(thresh_v_hig)
print 'Ending at points with V < ' + str(thresh_v_low)




#Now it gets the points of the regions
array_pos_region=array_results[(array_results[:,3] == 1),:]
#Resize of the box
l_x_max = np.max(array_pos_region[:,0])
l_y_max = np.max(array_pos_region[:,1])
l_z_max = np.max(array_pos_region[:,2])
l_x_min = np.min(array_pos_region[:,0])
l_y_min = np.min(array_pos_region[:,1])
l_z_min = np.min(array_pos_region[:,2])




fig = plt.figure(figsize=(25, 25))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(array_pos_region[:,0],array_pos_region[:,1],array_pos_region[:,2],s=20)
ax.set_xlim3d(l_x_min, l_x_max)
ax.set_ylim3d(l_y_min,l_y_max)
ax.set_zlim3d(l_z_min,l_z_max)
ax.set_xlabel('X kpc/h', size=30)
ax.set_ylabel('Y kpc/h', size=30)
ax.set_zlabel('Z kpc/h', size=30)
ax.set_title('Regions obtained ', size=50)
ax.tick_params(axis='both', which='major', labelsize=30)
plt.savefig("./regions_.png",format = 'png')






