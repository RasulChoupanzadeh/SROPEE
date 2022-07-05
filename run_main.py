""" run_main.py     => This is the main program of SROPEE. For implementing SROPEE. You only need to put the input variables in this script, and run it.
                                    

Author: Rasul Choupanzadeh 
Date: 07/05/2022


This program executes the  Full synthesis, Passive MOR, and Reduced synthesis parts of SROPEE, which are based on the concepts from [1-9].

 [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses
     by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

 [2] B. Gustavsen, "Improving the pole relocating properties of vector fitting",
     IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

 [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems
     Using a Fast Implementation of the Vector Fitting Method", IEEE Microwave and Wireless Components
     Letters, vol. 18, no. 6, pp. 383-385, June 2008.

 [4] A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization 
      through equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119.

 [5] V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from Frequency Response
      of Multi-port Networks", EMC EUR, pp. 79-84, 2016. [Online]. Available: ://WOS:000392194100012.

 [6] R. Choupanzadeh and A. Zadehgol. "Stability, causality, and passivity analysis of canonical equivalent 
      circuits of improper rational transfer functions with real poles and residues", IEEE Access, vol.8, pp. 125149-125162, 2020.
      
 [7] https://github.com/PhilReinhold/vectfit_python   

 [8] https://github.com/JenniferEHoule/Circuit_Synthesis
 
 [9] http://scikit-rf.org
 
 [10] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
      network analysis. IEEE Transactions on Circuits and Systems, 22(6), pp. 504-509, 1975.
    
 [11] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

 [12] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
      SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
      In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pp. 74-79, 2004.
    
 [13] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
      Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
      RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pp. 244-249, 2005.
    
"""



## Input: te_rough_example.s4p file (note: this file must be located in the folder: .\SROPEE\Input), Input_file_name, Num_pole_pairs, n, in_port            
## Output: input_variables.npy, full_netlist.sp, reduced_MNA_matrices.npz, reduced_netlist.sp, comparison figures for verification of obtained results in each part of SROPEE




import os
import numpy as np

#----------------------------------------------------------------Automated SROPEE inputs--------------------------------------------------------
Input_file_name = 'te_rough_example.s4p'

#Full synthesis input
Num_pole_pairs = 100

#Passive MOR input
n = 200                         # n = order of reduction, may be estimated according to the order of original network (N) => (e.g., n=0.4*N)
                                # N=2*Num_pole_pairs * number of sub-circuits, where Num_subckt=int((p*(p+1))/2)

#Reduced synthesis input
in_port = 1                    # port of input source


#Save input variables
np.save('./Output/input_variables.npy', [Input_file_name, Num_pole_pairs, n, in_port])

#--------------------------------------------------------------Executing Automated SROPEE algorithm------------------------------------------------------
import sys
sys.path.insert(0, '1_Full synthesis') 
sys.path.insert(1, '2_Passive MOR') 
sys.path.insert(2, '3_Reduced synthesis') 
import Run_full_synthesis
import Run_passive_MOR
import Run_reduced_synthesis


print('\n********************************* SROPEE algorithm is completed *************************************\n')

