
""" Run_main_passive_MOR.py      => This script is the main, which runs the entire passive model order reduction (passive MOR) and reduced synthesis.

Author: Rasul Choupanzadeh
Date: 06/28/2022


This code uses general_MNA_builder.py and Block_SAPOR.py, which are based on the concepts from [1-2] and [3-4], respectively.

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6), pp. 504-509, 1975.
    
[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

[3] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pp. 74-79, 2004.
    
[4] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pp. 244-249, 2005.
    
"""



## Input: Coef, fmin, fmax, fstep, n, s0, full_netlist.sp, Load_Z_ADS_full.py          Output: reduced_MNA_matrices.npz



import numpy as np
import matplotlib.pyplot as plt
import time
from general_MNA_builder import MNA_calculator
#from Load_equiv_circuit_results import Load_equiv_circuit_results
from sklearn.metrics import mean_squared_error
from math import sqrt

def precise_arange(start, stop, step):                             # Defining float range function with precise values
    return step * np.arange(start / step, stop / step)

        
#----------------------------------------------------------------------Input parameters-----------------------------------------------------------------------------------
Coef = 1e12                         # normalized coefficient in THz
fmin = 175
fmax = 215
fstep = 0.01
freq = precise_arange(fmin, fmax+fstep, fstep) 

n = 200                              # order of reduction
s0 = 2*np.pi*194*Coef

#load the netlist
exec(open("Load_netlist.py").read()) 

str_in = input('Please enter the desired Z-parameter in the form " Zmn " for m,n = 1,2,3,... (e.g., Z11, Z12,...):\n \n \t Zmn = ')
print()
in_port = int(str_in[2])-1
out_port = int(str_in[1])-1

Jm = np.zeros(shape=(Num_port,1))               
Jm[in_port] = 1                                                       

#-----------------------------------------------------------------General MNA approach--------------------------------------------------------------------------------
# build general MNA
[Gm,Cm,Gamma,Lv,Bm,M,N,SEt] = MNA_calculator(Ra,Rb,L,C,Num_branch, Num_port);
p = Bm.shape[1];
q = Lv.shape[1];


# Load full order network results (Z-parameter) from ADS---------------
## This part is for verification of original system and reduced system results with the simulated results of ADS software, and is disabled in final version.
#exec(open("Load_Z_ADS_full.py").read())
#Z_ADS_full = Z_ADS[:,out_port,in_port]
#plt.semilogy(freq, np.absolute(Z_ADS_full), label='ADS Full-order', c="black", linestyle='solid')


# Run original full-order network from General MNA matrices------------
exec(open("original_system.py").read())
t_full = 0;
t_full = t_full + end-start;
V_ports_full = V_orig[0:p,:]
Z_MNA_full = V_ports_full[out_port,:]/Jm[in_port]
plt.semilogy(freq, np.absolute(Z_MNA_full), label='MNA Full-order', c="blue", linestyle='dashdot')


# Run Block SAPOR Algorithm for Model Order Reduction------------------
exec(open("Block_SAPOR.py").read())
t_reduced = 0;
t_reduced = t_reduced + end-start;
V_ports_reduced = V_reduced[0:p,:]
Z_MNA_reduced = V_ports_reduced[out_port,:]/Jm[in_port]     
plt.semilogy(freq, np.absolute(Z_MNA_reduced), label='MNA Reduced-order', c="red", linestyle=(0,(5, 5)))        # (0,(5, 5)) is equal to a specific type of dashed
             
                
# Plot-----------------------------------------------------------------
plt.legend(loc='upper left')
plt.xlim([fmin,fmax])
plt.xlabel('Frequency (THz)')
plt.ylabel('$|Z|$')
plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)


# Extra info
print("The execution time of full-order MNA system is", t_full, 'seconds')
print("The execution time of reduced-order MNA system is", t_reduced, 'seconds')

RMSE = sqrt(mean_squared_error(np.absolute(Z_MNA_reduced),  np.absolute(Z_MNA_full)))
print("RMSE of reduced and full order MNA = ", RMSE)  

#RMSE = sqrt(mean_squared_error(np.absolute(Z_MNA_full),  np.absolute(Z_ADS_full)))
#print("RMSE of full-order MNA and ADS full-order network = ", RMSE) 


np.savez('reduced_MNA_data.npz', name1=Cmr, name2=Gmr, name3=Gammar, name4=Bmr, name5=Q, name6=Jm, name7=in_port, name8=out_port, name9=freq)

plt.show()
