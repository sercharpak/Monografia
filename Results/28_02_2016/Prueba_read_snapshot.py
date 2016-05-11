
# coding: utf-8

# In[1]:

#This module is from https://bitbucket.org/rthompson/pygadgetreader
from pygadgetreader import *


# In[39]:

#%pylab inline
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np



# In[3]:


snap_file_sim='../../Outputs/28_02_2016/snapshot_005'



requested='dmcount'
n_dm=readheader(snap_file_sim,requested)

print n_dm



#array_vel=readsnap(snap_file_example,'vel','dm')
#array_vel=readsnap(snap_file_sim,'vel','dm')
array_pos=readsnap(snap_file_sim,'pos','dm')
#array_pos=readsnap(snap_file_example,'pos','dm')

"""


l_x = 500000
l_y = 500000
l_z = 10000

array_pos_2=array_pos[((array_pos[:,0]<l_x) & (array_pos[:,1]<l_y) & (array_pos[:,2]<l_z)),:]


fig = plt.figure()
plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=10)
plt.savefig("./pos_plano.png",format = 'png')
plt.close(fig)

"""

array_vel=readsnap(snap_file_sim,'vel','dm')




array_magnitud_vel=np.zeros(n_dm)
for i in range(0,n_dm):
    array_magnitud_vel[i]=np.linalg.norm(array_vel[i,:])
print "La magnitud de las velocidades ha sido calculada \n"

"""
#---------------------------------------------------------------
fig = plt.figure()
binwidth=10
plt.hist(array_magnitud_vel, bins=np.arange(min(array_magnitud_vel), max(array_magnitud_vel) + binwidth, binwidth))
#plt.xlabel('
plt.title('Speed Histogram', fontsize=20)
plt.savefig("./hist_vel.png",format = 'png')
plt.close(fig)
"""
#---------------------------------------------------------------


print "Pasamos a los scatter con colores. El color depende de la magnitud de la velocidad. \n"


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
plt.savefig("./scatter_magnitud_vel"+str(l_z)+".png",format = 'png')

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
    plt.savefig("./scatters_color/scatter_magnitud_vel"+str(l_z_new)+".png",format = 'png')
    print "Ya termino el scatter plot con z entre "+str(l_z_old) +" y "+ str(l_z_new)
    l_z_old=l_z_new
    plt.close(fig)
    del array_pos_2
    del array_magnitud_vel_m


#---------------------------------------------------------------

minimo=np.min(array_magnitud_vel)
desviacion = np.std(array_magnitud_vel)
umbral_v=minimo+desviacion/2.0

print "El umbral de la velocidad es: "+str(umbral_v)+"\n"

array_pos_3=array_pos[(array_magnitud_vel[:]<(umbral_v)),:]
array_magnitud_vel_2=array_magnitud_vel[(array_magnitud_vel[:]<(umbral_v))]

print "Ya hizo el calculo el minimo y la desviacion. Ya formo los arreglos de posicion y magnitud velocidad \n"

print "El nuevo arreglo posiciones tiene forma \n"
print array_pos_3.shape

print "El nuevo arreglo velocidades tiene forma \n"
print array_magnitud_vel_2.shape



#---------------------------------------------------------------


print "Pasamos a los scatter con colores con los arreglos recien calculados. El color depende de la magnitud de la velocidad. \n"


l_x = 500000
l_y = 500000
l_z = 100000

array_pos_2=array_pos_3[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z)),:]
array_magnitud_vel_m=array_magnitud_vel_2[((array_pos_3[:,0]<l_x) & (array_pos_3[:,1]<l_y) & (array_pos_3[:,2]<l_z))]

fig=plt.figure(figsize=(20,20))
plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=200, c=array_magnitud_vel_m)
plt.tick_params(axis='both', which='major', labelsize=50)
plt.grid()
plt.title('Speed lesser than '+ str(umbral_v)+' \n for z between '+str(0) + ' and '+str(l_z),fontsize=60)
plt.xlabel('X kpc/h', fontsize=50)
plt.ylabel('Y kpc/h', fontsize=50)
cbar=plt.colorbar()
cbar.set_label('Speed',size=50)
cbar.ax.tick_params(labelsize=50) 
#plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(umbral_v)+"_lz_"+str(l_z)+".png",format = 'png')
plt.savefig("./scatter_vel_menor/scatter_magnitud_vel_"+str(umbral_v)+"_lz_"+str(l_z)+".png",format = 'png')

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
        plt.scatter(array_pos_2[:,0], array_pos_2[:,1],s=200, c=array_magnitud_vel_m)
        plt.tick_params(axis='both', which='major', labelsize=40)
        plt.grid()
        plt.title('Speed lesser than '+ str(umbral_v)+' \n for z between '+str(l_z_old) + ' and '+str(l_z),fontsize=60)
        plt.xlabel('X kpc/h', fontsize=50)
        plt.ylabel('Y kpc/h', fontsize=50)
        cbar=plt.colorbar()
        cbar.set_label('Speed',size=50)
        cbar.ax.tick_params(labelsize=40) 
        #plt.savefig("./scatter_vel_menor_500/scatter_magnitud_vel_"+str(umbral_v)+"_lz_"+str(l_z_new)+".png",format = 'png')
        plt.savefig("./scatter_vel_menor/scatter_magnitud_vel_"+str(umbral_v)+"_lz_"+str(l_z_new)+".png",format = 'png')
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
ax.set_xlabel('X kpc/h', size=50)
ax.set_ylabel('Y kpc/h', size=50)
ax.set_zlabel('Z kpc/h', size=50)
ax.set_title('Speed lesser than '+ str(umbral_v), size=60)
ax.tick_params(axis='both', which='major', labelsize=50)
#plt.colorbar()
#plt.show()
#plt.savefig("./pos_3d_vel_menor_"+str(minimo+desviacion/2.0)+".png",format = 'png')
plt.savefig("./pos_3d_vel_menor_"+str(umbral_v)+"_s_smaller.png",format = 'png')
plt.close(fig)

print "Ploteo en 3D a "+str(len(array_magnitud_vel_2))+" galaxias"

#---------------------------------------------------------------

maximo=np.max(array_magnitud_vel)
#desviacion = np.std(array_magnitud_vel)
umbral_v_max=maximo-desviacion/2.0

print "El umbral max de la velocidad es: "+str(umbral_v_max)+"\n"

array_pos_3=array_pos[(array_magnitud_vel[:]>(umbral_v_max)),:]
array_magnitud_vel_2=array_magnitud_vel[(array_magnitud_vel[:]>(umbral_v_max))]

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(array_pos_3[:,0],array_pos_3[:,1],array_pos_3[:,2],c=array_magnitud_vel_2,s=0.1)
ax.set_xlabel('X kpc/h', size=50)
ax.set_ylabel('Y kpc/h', size=50)
ax.set_zlabel('Z kpc/h', size=50)
ax.set_title('Speed greater than '+ str(umbral_v_max), size=60)
ax.tick_params(axis='both', which='major', labelsize=50)
#plt.colorbar()
#plt.show()
#plt.savefig("./pos_3d_vel_menor_"+str(minimo+desviacion/2.0)+".png",format = 'png')
plt.savefig("./pos_3d_vel_mayor_"+str(umbral_v_max)+"_s_smaller.png",format = 'png')
plt.close(fig)

print "Ploteo en 3D a "+str(len(array_magnitud_vel_2))+" galaxias"

#---------------------------------------------------------------

"""
print "Ya hizo el grafico general. Pasamos a los graficos de regiones."

l_x_ini = 0
l_y_ini = 0
l_z_ini = 0

l_x_old = l_x_ini
l_y_old = l_y_ini
l_z_old = l_z_ini

l_x_step = 250000
l_y_step = 250000
l_z_step = 250000

l_x_new = l_x_old + l_x_step
l_y_new = l_y_old + l_y_step
l_z_new = l_z_old + l_z_step

l_x_max = 500000
l_y_max = 500000
l_z_max = 500000

n_steps_x = (l_x_max/l_x_new)
n_steps_y = (l_y_max/l_y_new)
n_steps_z = (l_z_max/l_z_new)

for i in range(1,n_steps_x+1):
    for j in range(1,n_steps_y+1):
        for k in range(1,n_steps_z+1):
            array_pos_r=array_pos_3[((array_pos_3[:,0]<=l_x_new) & (array_pos_3[:,0]>=l_x_old) & (array_pos_3[:,1]<=l_y_new) & (array_pos_3[:,1]>=l_y_old) & (array_pos_3[:,2]<=l_z_new) & (array_pos_3[:,2]>=l_z_old)),:]
            array_magnitud_vel_r=array_magnitud_vel_2[((array_pos_3[:,0]<=l_x_new) & (array_pos_3[:,0]>=l_x_old) & (array_pos_3[:,1]<=l_y_new) & (array_pos_3[:,1]>=l_y_old) & (array_pos_3[:,2]<=l_z_new) & (array_pos_3[:,2]>=l_z_old))]
    
            print "El nuevo arreglo posiciones tiene forma \n"
            print array_pos_r.shape

            print "El nuevo arreglo velocidades tiene forma \n"
            print array_magnitud_vel_r.shape

            n_galaxias=len(array_magnitud_vel_r)
            if(n_galaxias>0):
                fig = plt.figure(figsize=(25, 25))
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(array_pos_r[:,0],array_pos_r[:,1],array_pos_r[:,2],c=array_magnitud_vel_r,s=10)
                ax.set_xlabel('X kpc/h', size=30)
                ax.set_ylabel('Y kpc/h', size=30)
                ax.set_zlabel('Z kpc/h', size=30)
                ax.set_title('Galaxias con |V| < '+str(minimo+desviacion/2.0)+' con: \n'+str(l_x_old)+'<x<'+str(l_x_new)+'\n'+str(l_y_old)+'<y<'+str(l_y_new)+'\n'+str(l_z_old)+'<z<'+str(l_z_new), size=30)
                ax.tick_params(axis='both', which='major', labelsize=50)
                plt.savefig("./scatter_3D_half/pos_3d_vel_menor_"+str(minimo+desviacion/2.0)+'_'+str(l_x_old)+'_x_'+str(l_x_new)+'_'+str(l_y_old)+'_y_'+str(l_y_new)+'_'+str(l_z_old)+'_z_'+str(l_z_new)+".png",format = 'png')
                plt.close(fig)
                print "Ploteo en 3D a "+str(n_galaxias)+" galaxias"
            l_z_old = l_z_new
            l_z_new = l_z_old + l_z_step
        l_z_old = l_z_ini
        l_z_new = l_z_old + l_z_step
        l_y_old = l_y_new
        l_y_new = l_y_old + l_y_step
    l_y_old = l_y_ini
    l_y_new = l_y_old + l_y_step
    l_x_old = l_x_new
    l_x_new = l_x_old + l_x_step
"""    
       



