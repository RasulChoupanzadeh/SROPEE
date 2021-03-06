
""" Run_full_synthesis.py     => This program is the main program of full synthesis part of SROPEE. In sum, it:
                                    - Creates a network (with the help of scikit-rf) for a given S-parameter
                                    - Enforces the reciprocity of created network by symmetrizing S-parameter matrix
                                    - Fits the Y-parameters of reciprocated network using vector fitting algorithm
                                    - Genrates a fitted network for comparison of fitted results with the results of original network
                                    - Generates a netlist representing equivalent circuit of full-order network 
                                    

Author: Rasul Choupanzadeh 
Date: 07/03/2022


This program uses vectfit.py and create_netlist.py program files, which are based on the concepts from [1-9].

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


"""


## Input: Input_file_name and Num_pole_pairs from input_variables.npy           Output: poles, residues, "full_netlist.sp" file, RMSE of fitted and actual network, comparison figures (network parameters) of actual and fitted network

import os

import vectfit
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from math import sqrt
from create_netlist import create_netlist_file

import skrf as rf
from skrf import Network, Frequency


def scattering_from_admittance(Y, z0):
    p = Y.shape[1]
    N_freq = Y.shape[0]
    S_par = np.zeros(shape=(N_freq, p, p), dtype='complex')
    for f in range( N_freq):
        E = np.eye(p)
        Zref = z0*E
        Gn = 1/sqrt(np.real(z0))
        Gref = Gn*E
        S_par[f,:,:] = Gref @ (E - Zref @ Y[f,:,:]) @ np.linalg.inv(E + Zref @ Y[f,:,:]) @ np.linalg.inv(Gref)

    return S_par

def admittance_from_scattering(S, z0):
    p = S.shape[1]
    N_freq = S.shape[0]
    Y_par = np.zeros(shape=(N_freq, p, p), dtype='complex')
    for f in range(N_freq):
        E = np.eye(p)
        Zref = z0*E
        Gn = 1/sqrt(np.real(z0))
        Gref = Gn*E
        Y_par[f,:,:] = np.linalg.inv(Gref) @ np.linalg.inv(S[f,:,:] @ Zref + Zref) @ (E - S[f,:,:]) @ Gref

    return Y_par


def create_upper_matrix(values, size):
    upper = np.zeros((size, size),dtype='complex')
    upper[np.triu_indices(size, 0)] = values
    return(upper)



# Where to save the figures
PROJECT_ROOT_DIR = "./Output/"
CHAPTER_ID = "Full_synthesis"
IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "figures", CHAPTER_ID)
os.makedirs(IMAGES_PATH, exist_ok=True)


def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)
    

#---------------------------------------------------------Load the origianl network results and enforce reciprocity--------------------------------------------------- 
Input_file_name = np.load('./Output/input_variables.npy')[0]
Num_pole_pairs = int(np.load('./Output/input_variables.npy')[1])

print('\n********************************* Starting part one (Full synthesis) *********************************\n')

input_full_synthesis = './Input/'+ Input_file_name
TL = Network(input_full_synthesis)
S_par_original = TL.s.copy()

# creating frequency matrix    
freq_start = TL.frequency.start   
freq_stop = TL.frequency.stop
Num_freq = TL.frequency.npoints
freq = np.linspace(freq_start,freq_stop,Num_freq)

number_of_ports = S_par_original.shape[1]


S_par_reciprocal = TL.s.copy()
# Enforce reciprocity 
for f in range(len(freq)):
    for i in range(number_of_ports):
        for j in range(number_of_ports):
            if j<i:
                S_par_reciprocal[f,i,j] = S_par_reciprocal[f,j,i]


Y_par_reciprocal = admittance_from_scattering(S_par_reciprocal, 50)


#---------------------------------------------------Prepare the reciprocal network (i.e., construct and stack y') for vectfit algorithm --------------------------------------- 
Y_par = Y_par_reciprocal

# Calculate Y', reshape and stack the Y-parameters
p = number_of_ports
Num_subckt = int((p*(p+1))/2)                                           # sum of first p natural numbers
Num_freq = len(freq)
Y_stack = np.zeros(shape=(Num_subckt, Num_freq), dtype='complex')
for f in range(len(freq)):
    t = Y_par[f,:,:].copy()
    # construct y' (for futher details see [5]) and stack 
    for i in range(p):
        t[i,i] = np.sum(t[i,:])
                
    for i in range(p):
        for j in range(p):
            if i != j:
                t[i,j] = -t[i,j]

    Y_stack[:,f] = t[np.triu_indices(p)]                                # converts upper triangular values into a vector 


#---------------------------------------------------------------------Perform vectfit algorithm----------------------------------------------------------------------- 
test_s = 1j*2*np.pi*freq
test_f = Y_stack.copy()

# Find the fitted poles and residues
poles = np.zeros(shape=(Num_subckt,2*Num_pole_pairs), dtype='complex')    
residues = np.zeros(shape=(Num_subckt,2*Num_pole_pairs), dtype='complex')    
for i in range(Num_subckt):
    poles[i,:], residues[i,:] = vectfit.vectfit_auto_rescale(test_f[i,:], test_s, n_poles=Num_pole_pairs)             # Note: n_poles = each pair of complex-conjugate poles are counted as one 


fitted = np.zeros(shape=(Num_subckt,len(test_s)), dtype='complex')    
for i in range(0,Num_subckt):
    fitted[i,:] = vectfit.model(test_s, poles[i,:], residues[i,:])

   
RMSE = []
for i in range(0,Num_subckt):
    square = np.square(np.absolute(fitted[i,:]-test_f[i,:]))
    MSE=square.mean()
    rmse = np.sqrt(MSE)
    RMSE.append(rmse)
    
print("Maximum RMSE of fitted (Python) and actual values = ", np.max(RMSE))


#------------------------------------------------Genrate fitted network and compare the results with original network (reciprocated)----------------------------------------- 
# Genrate fitted network 
Y_fit = np.zeros(shape=(len(freq),p,p), dtype='complex')
for f in range(len(freq)):
    u = create_upper_matrix(fitted[:,f], p)
    l = np.triu(u,1).T
    Y_fit[f,:,:] = u + l
    # convert from y' to y
    for i in range(p):
        Y_fit[f,i,i] = np.sum(Y_fit[f,i,:])
                
    for i in range(p):
        for j in range(p):
            if i != j:
                Y_fit[f,i,j] = -Y_fit[f,i,j]

S_par_fitted = scattering_from_admittance(Y_fit, 50)

Z_par_reciprocal = np.linalg.inv(Y_par_reciprocal)
Z_fit = np.linalg.inv(Y_fit)

def compare_fitted_vs_actual(m, n):
    mn = m*10+n
    # Compare Z-parameters
    actual_Z =  Z_par_reciprocal[:,m-1,n-1]
    fitted_Z = Z_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.plot(freq/1e12, np.absolute(actual_Z), color='blue', label='Actual', linestyle='solid')
    ax.plot(freq/1e12, np.absolute(fitted_Z), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Z)*180/np.pi, color='blue', label='Actual', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_Z)*180/np.pi, color='red', label='Fitted', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"Z%d" %mn)
    save_fig(f"Z%d" %mn)

    
    # Compare Y-parameters
    actual_Y =  Y_par_reciprocal[:,m-1,n-1]
    fitted_Y = Y_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.plot(freq/1e12, np.absolute(actual_Y), color='blue', label='Actual', linestyle='solid')
    ax.plot(freq/1e12, np.absolute(fitted_Y), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Y)*180/np.pi, color='blue', label='Actual', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_Y)*180/np.pi, color='red', label='Fitted', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"Y%d" %mn)
    save_fig(f"Y%d" %mn)


    # Compare S-parameters
    actual_S =  S_par_reciprocal[:,m-1,n-1]
    fitted_S = S_par_fitted[:,m-1,n-1]
    
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.plot(freq/1e12, np.absolute(actual_S), color='blue', label='Actual', linestyle='solid')
    ax.plot(freq/1e12, np.absolute(fitted_S), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_S)*180/np.pi, color='blue', label='Actual', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_S)*180/np.pi, color='red', label='Fitted', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"S%d" %mn)
    save_fig(f"S%d" %mn)


# Comparison plots 
for m in range(1,p+1):
    for n in range(m,p+1):
        compare_fitted_vs_actual(m,n)
        plt.close('all')
      


#-----------------------------------------------------------------------Genrate netlist------------------------------------------------------------------ 
create_netlist_file(poles, residues, number_of_ports)


print('\n********************************* Part one (Full synthesis) is done *********************************\n')
