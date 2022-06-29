
""" Load_Z_ADS_full.py      => This script Loads the Z-parameters of original full-order network calculated using Keysight ADS software

Author: Rasul Choupanzadeh
Date: 06/28/2022

"""

 
 ## Input: Z_param_full.s4p           Output: Z_ADS matrix
 
 
import numpy as np

# Load and read the text file line by line 
nstr = str(p)
file_name = 'Z_param_full' + '.s' + nstr + 'p'      
ADS = open(file_name).readlines()  

# Find the starting and ending rows of values
for i in range(len(ADS)):
    if ADS[i][0] == '!':
        start_value = i+1  
end_value = len(ADS)

ADS = ADS[start_value:end_value][:]                 # removes the initial text lines
ADS = [line.split() for line in ADS]                # removes the empty columns (space between numbers)
ADS = [ele for ele in ADS if ele != []]             # removes the empty rows (empty lines)

# remove the frequency values
for i in range(len(ADS)):
    if len(ADS[i]) == 2*p+1:
        ADS[i] = ADS[i][1:]

# converting from string to float
ADS_shape = np.shape(ADS)
Z_ADS_raw = np.zeros(shape=(ADS_shape[0],ADS_shape[1]))     
for i in range(len(ADS)):
    for j in range(len(ADS[i])):
        Z_ADS_raw[i][j] = float(ADS[i][j])

# Combining real and imaginary values into a single value
Z_ADS_raw2 = np.zeros(shape=(ADS_shape[0],p), dtype='complex')   
nt = 0
for k in range(0,2*p,2):
    Z_ADS_raw2[:,nt] = Z_ADS_raw[:,k]+1j*Z_ADS_raw[:,k+1]
    nt = nt+1

num_freq_points = int(Z_ADS_raw.shape[0]/p)
Z_ADS = np.zeros(shape=(num_freq_points,p,p), dtype='complex')
for i in range(0,num_freq_points):
    Z_ADS[i,:,:] = Z_ADS_raw2[i*p:(i+1)*p,:]

