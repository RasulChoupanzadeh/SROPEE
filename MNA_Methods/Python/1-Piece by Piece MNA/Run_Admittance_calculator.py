
''' Run_Admittance_calculator.py

Author: Rasul Choupanzadeh
Date: 05/12/2022


This code uses MNA_calculator.py and Block_SAPOR.py, which are based on the concepts from [1-2] and [3-4], respectively.

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

[3] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.

[4] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

'''



import numpy as np
import matplotlib.pyplot as plt
import time
from MNA_calculator import MNA_calculator
from Load_netlist import Load_sbck
from sklearn.metrics import mean_squared_error
from math import sqrt


def precise_arange(start, stop, step):                             # Defining float range function with precise values
    return step * np.arange(start / step, stop / step)

        
#----------------------------------------------------------------------Input parameters-----------------------------------------------------------------------------------
Coef = 1e12                         # normalized coefficient in THz
fmin = 175
fmax = 215
fstep = 0.1
freq = precise_arange(fmin, fmax+fstep, fstep)  
Jm = [1]

#load the netlist
exec(open("Load_netlist.py").read())

## fixed reduction order  ==> Instead of two following lines, you can define non-identical order of reduction; like n_arr=[n0,n1,...,n_sbck]
n0 = 35
n_arr = np.ones(Num_sbck)*n0        ## defines a fixed order of reduction for all sub-circuits n_arr=[n0, n0,...,n0]

s0 = 2*np.pi*194*Coef


str_in = input('Please enter the desired Y parameter in the form " Ymn ".\n \n \t Ymn = ')
print()
row = int(str_in[1])-1
col = int(str_in[2])-1

#-----------------------------------------------------------------Piece-by-Piece MNA approach--------------------------------------------------------------------------------
t_full = 0
t_reduced = 0
I_orig = np.zeros(shape=(Num_sbck,len(freq)), dtype='complex')
I_red = np.zeros(shape=(Num_sbck,len(freq)), dtype='complex')
for ks in range(Num_sbck):
    sbck = ks+1
    n = int(n_arr[ks])
    [Num_branch,Ra,L,Rb,C] = Load_sbck(sbck)
    [Gm,Cm,Gamma,Lv,Bm,M,N] = MNA_calculator(Ra,Rb,L,C,Num_branch)
    p = Bm.shape[1]
    q = Lv.shape[1]  
    exec(open("original_system.py").read())
    t_full = t_full + end-start
    I_sbck_k_orig = np.zeros(shape=(1,len(freq)), dtype='complex')
    for i in range(0, len(Ra)):
        I_sbck_k_orig = I_sbck_k_orig+ (V_orig[i,:]-Jm[0])/Ra[i]
    I_orig [ks,:]= I_sbck_k_orig[0,:]
    
    exec(open("Block_SAPOR.py").read())
    t_reduced = t_reduced + end-start
    
    I_sbck_k_red = np.zeros(shape=(1,len(freq)), dtype='complex')
    for i in range(0, len(Ra)):
        I_sbck_k_red = I_sbck_k_red+ (V_reduced[i,:]-Jm[0])/Ra[i]
    I_red [ks,:]= I_sbck_k_red[0,:]    
  

# Y-parameter calculation
Y_orig = np.zeros(shape=(len(freq),Num_port,Num_port), dtype='complex')
## upper triangular Y-parameters
nt = 0
for i in range(Num_port):
    for j in range(i,Num_port):
        Y_orig[:,i,j] = I_orig[nt,:]
        nt = nt+1

## lower triangular Y-parameters
for f in range(len(freq)):
    Y_orig[f,:,:] = Y_orig[f,:,:] + np.transpose(Y_orig[f,:,:]) - np.diag(np.diag(Y_orig[f,:,:]))
   
## diagonal Y-parameters 
for f in range(len(freq)):
    for i in range(Num_port):
        for j in range(0,Num_port):
            if i != j:
                Y_orig[f,i,i] = Y_orig[f,i,i]+Y_orig[f,i,j]
                
                
Y_reduced = np.zeros(shape=(len(freq),Num_port,Num_port), dtype='complex')
## upper triangular Y-parameters
nt = 0
for i in range(Num_port):
    for j in range(i,Num_port):
        Y_reduced[:,i,j] = I_red[nt,:]
        nt = nt+1

## lower triangular Y-parameters
for f in range(len(freq)):
    Y_reduced[f,:,:] = Y_reduced[f,:,:] + np.transpose(Y_reduced[f,:,:]) - np.diag(np.diag(Y_reduced[f,:,:]))
   
## diagonal Y-parameters 
for f in range(len(freq)):
    for i in range(Num_port):
        for j in range(0,Num_port):
            if i != j:
                Y_reduced[f,i,i] = Y_reduced[f,i,i]+Y_reduced[f,i,j]

            
# -------------------------------------------------------------------Select plot--------------------------------------------------------------------
# Plot
plt.semilogy(freq, np.absolute(Y_orig[:,row,col]), label='Full-order', c="black")
plt.semilogy(freq, np.absolute(Y_reduced[:,row,col]), label='Reduced-order', c="red", linestyle='dotted')

plt.legend(loc='upper left')
plt.xlim([fmin,fmax])
plt.xlabel('Frequency (THz)')
plt.ylabel('$|Y|$')
plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)

RMSE = sqrt(mean_squared_error(np.absolute(Y_orig[:,row,col]),  np.absolute(Y_reduced[:,row,col])))
print("The execution time of full-order system is :", t_full)
print("The execution time of reduced-order system is :", t_reduced)

print("RMSE=", RMSE)

plt.show()