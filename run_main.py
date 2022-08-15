
import os
import numpy as np


PROJECT_ROOT_DIR = "./Output/"
Output_PATH = os.path.join(PROJECT_ROOT_DIR)
os.makedirs(Output_PATH, exist_ok=True)

#----------------------------------------------------------------Automated SROPEE inputs--------------------------------------------------------
Input_file_name = 'THz_S_param.s4p'

#Full synthesis input
Num_poles = 50

#Passive MOR input
n = 144                         # n = order of reduction, may be estimated according to the order of original network (N) => (e.g., n=0.4*N)
                                # N=2*Num_pole_pairs * number of sub-circuits, where Num_subckt=int((p*(p+1))/2)

#Reduced synthesis input
in_port = 1                    # port of input source 1,2,...,p


#Save input variables
path = os.path.join(Output_PATH, 'input_variables.npy')
np.save(path, [Input_file_name, Num_poles, n, in_port])

#--------------------------------------------------------------Run Automated SROPEE algorithm------------------------------------------------------
import sys
sys.path.insert(0, '1_Full synthesis') 
sys.path.insert(1, '2_Passive MOR') 
sys.path.insert(2, '3_Reduced synthesis') 
import Run_full_synthesis
import Run_passive_MOR
import Run_reduced_synthesis


print('\n********************************* SROPEE algorithm is completed *************************************\n')
