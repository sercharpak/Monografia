#!/usr/bin/python

#Plotting the general properties of the detected groups
#By Sergio Daniel Hernandez Charpak

#-------------------------------------------
#Imports
import numpy as np
import matplotlib.pyplot as plt
import sys
#-------------------------------------------
# USAGE = "python pipeline_histograms.py file_groups"
USAGE = "python pipeline_histograms.py file_groups"

#-------------------------------------------
#Functions
#Function to extract the thresholds from a group filename
def get_thresholds(filename):
    file_Dat = filename.strip(".dat")
    file_array = file_Dat.split("_")
    thresh_search_FA = float(file_array[9])
    thresh_search_Trace = float(file_array[11])
    thresh_seeds_FA = float(file_array[4])
    thresh_seeds_Trace = float(file_array[6])
    return thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace

#Function to analyze the Inertia Tensor
#Based in http://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def inertia_tensor_case (I_a, I_b, I_c):
    case = 0
    #Case 1 Ia=Ib=Ic: spherical top;
    if(I_a==I_b==I_c):
        case = 1
        return case
    #Case 2 Ia=Ib<Ic: oblate symmetric top;
    elif(I_a==I_b < I_c):
        case = 2
        return case
    #Case 3 Ia<Ib=Ic: prolate symmetric top;
    elif(I_a < I_b == I_c):
        case = 3
        return case
    #Case 4 Ia<Ib<Ic: asymmetric top
    elif(I_a < I_b < I_c):
        case = 4
        return case
    #Case 0 - ?
    else:
        case = 0
        return case
#-------------------------------------------
if(len(sys.argv)!=2):
    print "Please use correctly"
    print USAGE
    sys.exit()

inputfile = sys.argv[1]
print "Processing the file: " + inputfile
thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(inputfile)
group_data = np.loadtxt(inputfile)
n_groups = len(group_data)
print "There are", n_groups,"groups"

#-------------------------------------------
#Volume distribution
volumes = group_data[:,7]

print "Volume distribution for \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace)
fig = plt.figure(figsize = (8,8))
binwidth=0.5
plt.hist(np.log10(volumes), bins=np.arange(min(np.log10(volumes)), max(np.log10(volumes)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume (number of cells) ")
plt.ylabel("Number of Groups")
plt.grid()
plt.title("Distribution of volume (number of cells) with "+str(n_groups) + " groups \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace))
plt.savefig("./volumes_distr_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#-------------------------------------------
#Inertia Tensor Analysis
#Constants for writting the result of the shape of the respective detected group
CASES = []
CASE_0 = "none of the above"
CASES.append(CASE_0)
CASE_1 = "spherical top"
CASES.append(CASE_1)
CASE_2 = "oblate symmetric top"
CASES.append(CASE_2)
CASE_3 = "prolate symmetric top"
CASES.append(CASE_3)
CASE_4 = "asymmetric top"
CASES.append(CASE_4)

lambdas_a = group_data[:,0]
lambdas_b = group_data[:,1]
lambdas_c = group_data[:,2]

#This counting could be done way better
inertia_cases = np.zeros(5)

cases = np.zeros(n_groups)
for i in range(n_groups):
    cases [i] = inertia_tensor_case(lambdas_a[i], lambdas_b[i], lambdas_c[i])
    inertia_cases[cases [i]] += 1

print "Inertia Analysis for \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace)
fig = plt.figure(figsize = (8,8))
plt.hist(cases)
plt.xlabel("Inertia Case")
plt.ylabel("Number of Groups")
plt.grid()
plt.xlim(0,4.1)
plt.xticks([0.0,1.0,2.0,3.0,4.0])
plt.title("Inertia Cases for groups with "+str(n_groups) + " groups \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace))
plt.savefig("./inertia_cases_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

print "-------------------------------------------"
print "Case", "Case Number", "Number of occurences"
for i in range (len(inertia_cases)):
    print str(CASES[i]), str(i), str(inertia_cases[i])
