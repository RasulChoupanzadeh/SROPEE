'''
Author: Rasul Choupanzadeh
Date: 08/23/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

'''

import os
import numpy as np


PROJECT_ROOT_DIR = "./Output/"
Output_PATH = os.path.join(PROJECT_ROOT_DIR)
os.makedirs(Output_PATH, exist_ok=True)

#----------------------------------------------------------------Automated SROPEE inputs--------------------------------------------------------
################################################
# General input
################################################
Input_file_name = 'THz_S_param.s4p'


################################################
# Full synthesis inputs
################################################
Num_poles = 50

VF_iter1 = 4                            # Iteration number for calculating improved initial poles by fitting column sum in vector fitting algorithm
VF_iter2 = 10                           # Iteration number for calculating poles and residues in vector fitting algorithm
weight_f = 'sqrt_abs'                   # Choose between 'unweighted', 'abs', and 'sqrt_abs' weight types for f                      Note: 'abs' is strongest inverse wieght
weight_column_sum = 'norm'              # Choose between 'unweighted', 'norm', and 'sqrt_norm' weight types for g (column sum)       Note: 'norm' is strongest inverse wieght
VF_relax = True                         # Enable/disable relaxed vector fitting
VF_stable = True                        # Enable/disable enforce stablility of fitted poles
VF_asymp = 1                            # Fitting model (asymp=0  fits with None, asymp=1  fits with D, asymp=2  fits with D and E)

Passivity_Enforcement_Enable = True     # Enable/disable passivity assessment/enforcement 
SMP_iter_upper_limit = 10               # Passivity enforcement iteration upper limit

################################################
# Passive MOR input
################################################
n = 144                         # n = order of reduction, may be estimated according to the order of original network (N)
                                # N=2*Num_pols (i.e.,p) * number of sub-circuits, where Num_subckt=int((p*(p+1))/2)

################################################
# Reduced synthesis input
################################################
in_port = 1                     # port of input source 1,2,...,p




# Saving all input variables into .npy file
path = os.path.join(Output_PATH, 'input_variables.npy')
np.save(path, [Input_file_name, Num_poles, VF_iter1, VF_iter2, weight_f, weight_column_sum, VF_relax, VF_stable, VF_asymp, Passivity_Enforcement_Enable, SMP_iter_upper_limit, n, in_port])


#--------------------------------------------------------------Run Automated SROPEE algorithm------------------------------------------------------
import sys
sys.path.insert(0, '1_Full synthesis') 
sys.path.insert(1, '2_Passive MOR') 
sys.path.insert(2, '3_Reduced synthesis') 
import Run_full_synthesis
import Run_passive_MOR
import Run_reduced_synthesis


print('\n********************************* SROPEE algorithm is completed *************************************\n')
