#!/usr/bin/python
#
#Written by Sergio Hernandez Charpak
#

# Usage
# python pipeline_script_1.py file_eigen_vale_1 file_eigen_vale_2 file_eigen_vale_3 FA_thresh_value Div_thresh_value
#

USAGE = "python pipeline_script_1.py file_eigen_vale_1 file_eigen_vale_2 file_eigen_vale_3 FA_thresh_value Div_thresh_value"

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

#Function to read the EigenValues from the v_web CIC outpu
def read_CIC_scalar(filename):
    f = open(filename, "rb")
    dumb = f.read(38)

    dumb = f.read(4)
    n_x = f.read(4)
    n_y = f.read(4)
    n_z = f.read(4)
    nodes = f.read(8)
    x0 = f.read(4)
    y0 = f.read(4)
    z0 = f.read(4)
    dx = f.read(4)
    dy = f.read(4)
    dz = f.read(4)
    dumb = f.read(4)

    n_x = (unpack('i', n_x))[0]
    n_y = (unpack('i', n_y))[0]
    n_z = (unpack('i', n_z))[0]
    nodes = (unpack('q', nodes))[0]
    dx = (unpack('f', dx))[0]
    dy = (unpack('f', dy))[0]
    dz = (unpack('f', dz))[0]
    x0 = (unpack('f', x0))[0]
    y0 = (unpack('f', y0))[0]
    z0 = (unpack('f', z0))[0]
    print n_x, n_y, n_z, nodes, dx, dy, dz

    total_nodes = n_x * n_y *n_z
    dumb = f.read(4)
    array_data = f.read(total_nodes*4)
    dumb = f.read(4)
    format_s = str(total_nodes)+'f'
    array_data = unpack(format_s, array_data)
    f.close()
    array_data  = np.array(array_data)
    array_data.resize(n_z,n_y,n_x)
    array_data = array_data.transpose()
    return array_data

#---------------------------------------------------------

if(len(sys.argv)!=6):
    print "Please use correctly"
    print USAGE
    sys.exit() 

inputfile_1 = sys.argv[1]
inputfile_2 = sys.argv[2]
inputfile_3 = sys.argv[3]

print "Preparing to read the eigenvalues"
print inputfile_1,inputfile_2,inputfile_3

grid_1 = read_CIC_scalar(inputfile_1)
grid_2 = read_CIC_scalar(inputfile_2)
grid_3 = read_CIC_scalar(inputfile_3)

print "Eigenvalues read"

print "Forming the FA and the Trace grids"

FA = (grid_1-grid_3)**2  + (grid_2-grid_3)**2  + (grid_1-grid_2)**2 
FA = FA/(grid_1**2 + grid_2**2 + grid_3**2)
FA = np.sqrt(FA)/np.sqrt(3.0)

Trace = grid_1 + grid_2 + grid_3

print "Fa and Trace matrix formed"

thresh_FA = float(sys.argv[4])
thresh_Trace = float(sys.argv[5])

print "The Thresholds are: FA, Trace", str(thresh_FA), str(thresh_Trace) 

n_x,n_y,n_z = shape(FA)

print "Forming the resulting grid"

cell_list = []
cells_total = 0
for i in range (n_x):
    for j in range (n_y): 
        for k in range (n_z):
            if((FA[i,j,k]<thresh_FA ) & (Trace[i,j,k]>thresh_Trace)):
                cells_total +=1
                cell_list.append([i,j,k])

print "Grid formed with number of cells:", cells_total

output_name = '../FoF/src/grid_FA_'+str(thresh_FA)+'_Trace_'+str(thresh_Trace)+'.dat'

print "Writing in file " + output_name

fileout = open(output_name, 'w')

fileout.write("%d\n"%cells_total) #points in total 
fileout.write("%d\n"%cells_total) #points in 'DM'
fileout.write("0\n") #gas
fileout.write("0\n") #stars
fileout.write("0.01\n") # time
fileout.write("0\n") # nactive

for l in range(cells_total):
    fileout.write("%f %f %f\n"%(cell_list[l][0],cell_list[l][1],cell_list[l][2]))

fileout.close()

print "Finished writting the file"
