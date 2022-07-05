
""" Load_equiv_circuit_results.py      => This script Loads the results of equivalent reduced-order network simulated in Keysight ADS software which is exported as "reduced_voltages_ADS.cti" file

Author: Rasul Choupanzadeh
Date: 07/03/2022

"""

 
## Input: n, Q, and "reduced_voltages_ADS.cti" file           Output: V_equiv
 

import numpy as np

def Load_equiv_circuit_results(n,Q):
    
    # Load and read the netlist line by line 
    equiv_results = open('./Output/reduced_voltages_ADS.cti').readlines() 
    equiv_shape = np.shape(equiv_results)
    
    # Find the number of frequency points
    for line in equiv_results:
        if 'freq' in line:
            res = [int(i) for i in line.split() if i.isdigit()]                 # Extract integer numbers from string
            N_freq = int(res[0])
    
    V_r = np.zeros(shape=(n,N_freq), dtype= 'complex') 
    nr = 0
    for i in range(equiv_shape[0]):
        if  equiv_results[i] == 'BEGIN\n':
            start_value = i + 1
            stop_value = i + N_freq+1
            nc = 0
            for k in range(start_value,stop_value):
                res = equiv_results[k].split(',')
                V_r[nr,nc] = float(res[0]) + 1j*float(res[1])
                nc = nc + 1
            nr = nr + 1
        
    V_equiv = Q @ V_r
            
    return V_equiv

