from pygadgetreader import *
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
matplotlib.use('Agg') 

snap_file_sim='/lustre/home/ciencias/fisica/pregrado/sd.hernandez204/Outputs/28_02_2016/snapshot_005'
results_path = '/lustre/home/ciencias/fisica/pregrado/sd.hernandez204/Outputs/graphics/09_03_2016/'

scatter_general_path = results_path + 'scatter_general/'
scatter_low_path = results_path + 'scatter_low/'
scatter_high_path = results_path + 'scatter_high/'


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


fig = plt.figure(figsize=(20,20))
binwidth=10
plt.hist(array_magnitud_vel, bins=np.arange(min(array_magnitud_vel), max(array_magnitud_vel) + binwidth, binwidth))
#plt.xlabel('
plt.title('Speed Histogram', fontsize=20)
plt.savefig(results_path+"hist_vel.png",format = 'png')
plt.close(fig)


print "Pasamos a los scatter con colores. El color depende de la magnitud de la velocidad. \n"


vel_max = np.max(array_magnitud_vel)
vel_min = np.min(array_magnitud_vel)
vel_std = np.std(array_magnitud_vel)
vel_mean = np.mean(array_magnitud_vel)

print "Vel_max, Vel_min, std_vel, vel_mean"
print vel_max, vel_min, vel_std, vel_mean
thresh_v_high = vel_mean + 2.0 * vel_std
thresh_v_low = vel_mean -  vel_std
print "Thresh_high, Thresh_low"
print thresh_v_high, thresh_v_low


l_x = 500000
l_y = 500000
l_z = 100000

array_pos_2=array_pos[((array_pos[:,0]<l_x) & (array_pos[:,1]<l_y) & (array_pos[:,2]<l_z)),:]
array_magnitud_vel_2=array_magnitud_vel[((array_pos[:,0]<l_x) & (array_pos[:,1]<l_y) & (array_pos[:,2]<l_z))]

fig=plt.figure(figsize=(20,20))
plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_2)
plt.tick_params(axis='both', which='major', labelsize=40)
plt.grid()
plt.title('Speed distribution for \n  z in [ '+str(0) + ' ; '+str(l_z)+' ]',fontsize=60)
plt.xlabel('X kpc/h', fontsize=50)
plt.ylabel('Y kpc/h', fontsize=50)
cbar=plt.colorbar()
cbar.set_label('Speed',size=50)
cbar.ax.tick_params(labelsize=40) 
plt.savefig(scatter_general_path+"scatter_magnitud_vel"+str(l_z)+".png",format = 'png')

print "Ya termino el scatter plot con z entre "+str(0) +" y "+ str(l_z)
del array_pos_2
del array_magnitud_vel_2
del fig
l_z_old=l_z
for i in range(1,5):
    l_z_new=l_z+i*l_z
    array_pos_2=array_pos[((array_pos[:,0]<l_x) & (array_pos[:,1]<l_y) & (array_pos[:,2]<l_z_new) & (array_pos[:,2]>l_z_old)),:]
    array_magnitud_vel_m=array_magnitud_vel[((array_pos[:,0]<l_x) & (array_pos[:,1]<l_y) & (array_pos[:,2]<l_z_new) & (array_pos[:,2]>l_z_old))]
    fig=plt.figure(figsize=(20,20))
    plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_m)
    plt.tick_params(axis='both', which='major', labelsize=40)
    plt.grid()
    plt.title('Speed distribution for \n  z in [ '+str(l_z_old) + ' ; '+str(l_z)+' ]',fontsize=60)
    plt.xlabel('X kpc/h', fontsize=50)
    plt.ylabel('Y kpc/h', fontsize=50)
    cbar=plt.colorbar()
    cbar.set_label('Speed',size=50)
    cbar.ax.tick_params(labelsize=40) 
    plt.savefig(scatter_general_path+"scatter_magnitud_vel"+str(l_z_new)+".png",format = 'png')
    print "Ya termino el scatter plot con z entre "+str(l_z_old) +" y "+ str(l_z_new)
    l_z_old=l_z_new
    plt.close(fig)
    del array_pos_2
    del array_magnitud_vel_m


#---------------------------------------------------------------



print "El umbral de la velocidad es: "+str(thresh_v_low)

array_pos_3=array_pos[(array_magnitud_vel[:]<(thresh_v_low)),:]
array_magnitud_vel_2=array_magnitud_vel[(array_magnitud_vel[:]<(thresh_v_low))]

print "Ya hizo el calculo el minimo y la desviacion. Ya formo los arreglos de posicion y magnitud velocidad "

print "El nuevo arreglo posiciones tiene forma"
print array_pos_3.shape

print "El nuevo arreglo velocidades tiene forma"
print array_magnitud_vel_2.shape



#---------------------------------------------------------------


print "Pasamos a los scatter con colores con los arreglos recien calculados. El color depende de la magnitud de la velocidad. \n"


l_x = 500000
l_y = 500000
l_z = 100000

array_pos_2=array_pos_3[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z)),:]
array_magnitud_vel_m=array_magnitud_vel_2[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z))]

fig=plt.figure(figsize=(20,20))
plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_m)
plt.tick_params(axis='both', which='major', labelsize=35)
plt.grid()
plt.title('Speed lesser than '+ str(thresh_v_low)+' \n for z between '+str(0) + ' and '+str(l_z),fontsize=40)
plt.xlabel('X kpc/h', fontsize=40)
plt.ylabel('Y kpc/h', fontsize=40)
cbar=plt.colorbar()
cbar.set_label('Speed',size=40)
cbar.ax.tick_params(labelsize=40) 
#plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z)+".png",format = 'png')
plt.savefig(scatter_low_path+"scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z)+".png",format = 'png')

print "Ya termino el scatter plot con z entre "+str(0) +" y "+ str(l_z)
del array_pos_2
del array_magnitud_vel_m
plt.close(fig)
l_z_old=l_z
#for i in range(1,500):
for i in range(1,5):
    l_z_new=l_z+i*l_z
    array_pos_2=array_pos_3[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z_new) & (array_pos_3[:,2]>l_z_old)),:]
    array_magnitud_vel_m=array_magnitud_vel_2[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z_new) & (array_pos_3[:,2]>l_z_old))]
    print "El nuevo arreglo posiciones tiene forma \n"
    print array_pos_2.shape

    print "El nuevo arreglo velocidades tiene forma \n"
    print array_magnitud_vel_m.shape

    n_galaxies=len(array_magnitud_vel_m)
    if(n_galaxies>0):
        
        fig=plt.figure(figsize=(20,20))
        plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_m)
        plt.tick_params(axis='both', which='major', labelsize=35)
        plt.grid()
        plt.title('Speed lesser than '+ str(thresh_v_low)+' \n for z between '+str(l_z_old) + ' and '+str(l_z),fontsize=40)
        plt.xlabel('X kpc/h', fontsize=40)
        plt.ylabel('Y kpc/h', fontsize=40)
        cbar=plt.colorbar()
        cbar.set_label('Speed',size=40)
        cbar.ax.tick_params(labelsize=40) 
        #plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z_new)+".png",format = 'png')
        plt.savefig(scatter_low_path+"scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z_new)+".png",format = 'png')
        print "Ya termino el scatter plot con z entre "+str(l_z_old) +" y "+ str(l_z_new)
        plt.close(fig)
        del fig
    l_z_old=l_z_new
    
    del array_pos_2
    del array_magnitud_vel_m

#--------------------------------------------------------------

print "Pasa a los plots en 3D \n"

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(array_pos_3[:,0],array_pos_3[:,1],array_pos_3[:,2],c=array_magnitud_vel_2,s=0.1)
ax.set_xlabel('X kpc/h', size=40)
ax.set_ylabel('Y kpc/h', size=40)
ax.set_zlabel('Z kpc/h', size=40)
ax.set_title('Particles with speed lesser than '+ str(thresh_v_low), size=40)
ax.tick_params(axis='both', which='major', labelsize=40)
#plt.colorbar()
#plt.show()
#plt.savefig("./pos_3d_vel_menor_"+str(minimo+desviacion/2.0)+".png",format = 'png')
plt.savefig(results_path+"pos_3d_vel_menor_"+str(thresh_v_low)+"_s_smaller.png",format = 'png')
plt.close(fig)

print "Ploteo en 3D a "+str(len(array_magnitud_vel_2))+" galaxias"

#--------------------------------------------------------------

print "El umbral de la velocidad es: "+str(thresh_v_high)

array_pos_3=array_pos[(array_magnitud_vel[:]>(thresh_v_high)),:]
array_magnitud_vel_2=array_magnitud_vel[(array_magnitud_vel[:]>(thresh_v_high))]

print "Ya hizo el calculo el minimo y la desviacion. Ya formo los arreglos de posicion y magnitud velocidad"

print "El nuevo arreglo posiciones tiene forma"
print array_pos_3.shape

print "El nuevo arreglo velocidades tiene forma"
print array_magnitud_vel_2.shape

l_x = 500000
l_y = 500000
l_z = 100000

array_pos_2=array_pos_3[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z)),:]
array_magnitud_vel_m=array_magnitud_vel_2[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z))]

fig=plt.figure(figsize=(20,20))
plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_m)
plt.tick_params(axis='both', which='major', labelsize=35)
plt.grid()
plt.title('Speed greater than '+ str(thresh_v_high)+' \n for z between '+str(0) + ' and '+str(l_z),fontsize=40)
plt.xlabel('X kpc/h', fontsize=40)
plt.ylabel('Y kpc/h', fontsize=40)
cbar=plt.colorbar()
cbar.set_label('Speed',size=40)
cbar.ax.tick_params(labelsize=40) 
#plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z)+".png",format = 'png')
plt.savefig(scatter_high_path+"scatter_magnitud_vel_"+str(thresh_v_high)+"_lz_"+str(l_z)+".png",format = 'png')

print "Ya termino el scatter plot con z entre "+str(0) +" y "+ str(l_z)
del array_pos_2
del array_magnitud_vel_m
plt.close(fig)
l_z_old=l_z
#for i in range(1,500):
for i in range(1,5):
    l_z_new=l_z+i*l_z
    array_pos_2=array_pos_3[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z_new) & (array_pos_3[:,2]>l_z_old)),:]
    array_magnitud_vel_m=array_magnitud_vel_2[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z_new) & (array_pos_3[:,2]>l_z_old))]
    print "El nuevo arreglo posiciones tiene forma "
    print array_pos_2.shape

    print "El nuevo arreglo velocidades tiene forma "
    print array_magnitud_vel_m.shape

    n_galaxies=len(array_magnitud_vel_m)
    if(n_galaxies>0):
        
        fig=plt.figure(figsize=(20,20))
        plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=20, c=array_magnitud_vel_m)
        plt.tick_params(axis='both', which='major', labelsize=35)
        plt.grid()
        plt.title('Speed greater than '+ str(thresh_v_high)+' \n for z between '+str(l_z_old) + ' and '+str(l_z),fontsize=40)
        plt.xlabel('X kpc/h', fontsize=40)
        plt.ylabel('Y kpc/h', fontsize=40)
        cbar=plt.colorbar()
        cbar.set_label('Speed',size=40)
        cbar.ax.tick_params(labelsize=40) 
        #plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(thresh_v_low)+"_lz_"+str(l_z_new)+".png",format = 'png')
        plt.savefig(scatter_high_path+"scatter_magnitud_vel_"+str(thresh_v_high)+"_lz_"+str(l_z_new)+".png",format = 'png')
        print "Ya termino el scatter plot con z entre "+str(l_z_old) +" y "+ str(l_z_new)
        plt.close(fig)
        del fig
    l_z_old=l_z_new
    
    del array_pos_2
    del array_magnitud_vel_m

#---------------------------------------------------------------

print "Pasa a los plots en 3D \n"

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(array_pos_3[:,0],array_pos_3[:,1],array_pos_3[:,2],c=array_magnitud_vel_2,s=0.1)
ax.set_xlabel('X kpc/h', size=40)
ax.set_ylabel('Y kpc/h', size=40)
ax.set_zlabel('Z kpc/h', size=40)
ax.set_title('Particles with speed greater than '+ str(thresh_v_high), size=40)
ax.tick_params(axis='both', which='major', labelsize=40)
#plt.colorbar()
#plt.show()
#plt.savefig("./pos_3d_vel_menor_"+str(minimo+desviacion/2.0)+".png",format = 'png')
plt.savefig(results_path+"pos_3d_vel_mayor_"+str(thresh_v_high)+"_s_smaller.png",format = 'png')
plt.close(fig)

print "Ploteo en 3D a "+str(len(array_magnitud_vel_2))+" galaxias"

