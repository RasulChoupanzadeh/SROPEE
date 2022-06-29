
''' Run_Admittance_calculator.py

Author: Rasul Choupanzadeh
Date: 05/12/2022


This code uses general_MNA_builder.py and Block_SAPOR.py, which are based on the concepts from [1-2] and [3-4], respectively.

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
from general_MNA_builder import MNA_calculator
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

n = 200                              # order of reduction
s0 = 2*np.pi*194*Coef

#load the netlist
exec(open("Load_netlist.py").read())

str_in = input('Please enter the desired Y parameter in the form " Ymn ".\n \n \t Ymn = ')
print()
row = int(str_in[1])-1
col = int(str_in[2])-1

Jm = np.zeros(shape=(Num_port,1))               
Jm[col] = 1                                                       

#-----------------------------------------------------------------General MNA approach--------------------------------------------------------------------------------
# build general MNA
[Gm,Cm,Gamma,Lv,Bm,M,N,SEt] = MNA_calculator(Ra,Rb,L,C,Num_branch, Num_port);

# Run full-order network from General MNA matrices------------------
p = Bm.shape[1];
q = Lv.shape[1];
exec(open("original_system.py").read())
t_full = 0;
t_full = t_full + end-start;

V_node = V_orig[p:,:]
I_ind = np.zeros(shape=(M,len(freq)), dtype='complex')
nt = 0
for f in freq:
    s = 1j*2*np.pi*f*Coef
    I_ind[:,nt] = (SEt[:,p:] @ V_node[:,nt])/s
    nt = nt+1
    
I_sbck = np.zeros(shape=(Num_sbck,len(freq)), dtype='complex')
nt = 0
for i in range(0,Num_sbck):
    I_sbck[i,:] = np.sum(I_ind[nt:nt+Num_branch[i],:],axis=0)
    nt= nt + Num_branch[i]


Y_orig = np.zeros(shape=(len(freq),p,p), dtype='complex')
## upper triangular Y-parameters
nt = 0
for i in range(p):
    for j in range(i,p):
        Y_orig[:,i,j] = I_sbck[nt,:]
        nt = nt+1

## lower triangular Y-parameters
for f in range(len(freq)):
    Y_orig[f,:,:] = Y_orig[f,:,:] + np.transpose(Y_orig[f,:,:]) - np.diag(np.diag(Y_orig[f,:,:]))

            
## diagonal Y-parameters 
for f in range(len(freq)):
    for i in range(p):
        for j in range(i+1,p):
            Y_orig[f,i,i] = Y_orig[f,i,i]+Y_orig[f,i,j]
            
        for j in range(0,i):
            Y_orig[f,i,i] = Y_orig[f,i,i]-Y_orig[f,i,j]        
                


# Run reduced-order network from General MNA matrices---------------
exec(open("Block_SAPOR.py").read())
t_reduced = 0;
t_reduced = t_reduced + end-start;

V_node_reduced = V_reduced[p:,:]
I_ind_reduced = np.zeros(shape=(M,len(freq)), dtype='complex')
nt = 0
for f in freq:
    s = 1j*2*np.pi*f*Coef
    I_ind_reduced[:,nt] = (SEt[:,p:] @ V_node_reduced[:,nt])/s
    nt = nt+1
    
I_sbck_reduced = np.zeros(shape=(Num_sbck,len(freq)), dtype='complex')
nt = 0
for i in range(0,Num_sbck):
    I_sbck_reduced[i,:] = np.sum(I_ind_reduced[nt:nt+Num_branch[i],:],axis=0)
    nt= nt + Num_branch[i]    

Y_reduced = np.zeros(shape=(len(freq),p,p), dtype='complex')
## upper triangular Y-parameters
nt = 0
for i in range(p):
    for j in range(i,p):
        Y_reduced[:,i,j] = I_sbck_reduced[nt,:]
        nt = nt+1

## lower triangular Y-parameters
for f in range(len(freq)):
    Y_reduced[f,:,:] = Y_reduced[f,:,:] + np.transpose(Y_reduced[f,:,:]) - np.diag(np.diag(Y_reduced[f,:,:]))

            
## diagonal Y-parameters 
for f in range(len(freq)):
    for i in range(p):
        for j in range(i+1,p):
            Y_reduced[f,i,i] = Y_reduced[f,i,i]+Y_reduced[f,i,j]
            
        for j in range(0,i):
            Y_reduced[f,i,i] = Y_reduced[f,i,i]-Y_reduced[f,i,j]        
                

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
