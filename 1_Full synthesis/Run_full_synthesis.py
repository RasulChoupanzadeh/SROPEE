
""" Run_full_synthesis.py     => This program is the main program of Full Synthesis part of SROPEE. In sum, it:
                                    - Creates a network (with the help of scikit-rf) for a given S-parameter
                                    - Enforces the reciprocity of created network by symmetrizing S-parameter matrix
                                    - Fits the Y-parameters of reciprocated network using vector fitting algorithm
                                    - Genrates a fitted network for comparison of fitted results with the results of original network
                                    - Assess/Enforce the passivity using singular test matrix 
                                    - Generates a netlist representing equivalent circuit of full-order network 
                                    

Author: Rasul Choupanzadeh 
Date: 08/10/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [12].

This program uses vectfit3.py and create_netlist.py program files, which are based on the concepts from [1-11].

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
      
 [7] Simon De Ridder, GitHub. Feb 08, 2019. Accessed on: August 10, 2022, [Online].https://github.com/SimonDeRidder/pyVF  


 [8] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: August 11, 2022, [Online]. https://github.com/JenniferEHoule/Circuit_Synthesis


 [9] Totorica, Nathan, GitHub, May 5, 2021. Accessed on: August 10, 2022, [Online]. Available:https://github.com/ntotorica/SMP_Passivity_Enforcement
 
 
 [10] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery,
    vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

 
 [11] http://scikit-rf.org
 
 [12] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.


"""


## Input: Input_file_name, Num_poles, and options from input_variables.npy           Output: poles, residues, "full_netlist.sp" file, RMSE of fitted and actual network, comparison figures (network parameters) of actual and fitted network

import os

import numpy as np
nax = np.newaxis
from vectfit3 import vectfit3, tri2full, ss2pr
from time import time
import matplotlib.pyplot as plt
from pylab import *
from math import sqrt
from eig_plot import plot
from create_netlist import create_netlist_file
from SMP import SMP
import skrf as rf



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


def convert_y_to_y_prime(y):
    number_of_ports = y.shape[1]
    y_prime = np.zeros(shape=y.shape, dtype='complex')
    for port_1 in range(0, number_of_ports):
        for port_2 in range(0, number_of_ports):
            if port_1 == port_2:
                y_prime[:, port_1, port_2] = np.sum(y[:, port_1, :], axis=1)
            else:
                y_prime[:, port_1, port_2] = -y[:, port_1, port_2].copy()
    return y_prime


# Where to save the figures
PROJECT_ROOT_DIR = "./Output/"
CHAPTER_ID = "Full_synthesis"
SECTION_ID1 = "1_Vector Fitting"
SECTION_ID2 = "2_Passivity Enforcement"
IMAGES_PATH1 = os.path.join(PROJECT_ROOT_DIR, "figures", CHAPTER_ID,SECTION_ID1)
IMAGES_PATH2 = os.path.join(PROJECT_ROOT_DIR, "figures", CHAPTER_ID,SECTION_ID2)
os.makedirs(IMAGES_PATH1, exist_ok=True)
os.makedirs(IMAGES_PATH2, exist_ok=True)



def save_fig_VF(fig_id, tight_layout=True, fig_extension="pdf", resolution=300):
    path = os.path.join(IMAGES_PATH1, fig_id + "." + fig_extension)
    print("Saving VF figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)
    
    
def save_fig_SMP(fig_id, tight_layout=True, fig_extension="pdf", resolution=300):
    path = os.path.join(IMAGES_PATH2, fig_id + "." + fig_extension)
    print("Saving Passive VF figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)    

#---------------------------------------------------------Load the origianl network results and enforce reciprocity--------------------------------------------------- 
Input_file_name = np.load('./Output/input_variables.npy')[0]
Num_poles = int(np.load('./Output/input_variables.npy')[1])
VF_iter1 = int(np.load('./Output/input_variables.npy')[2])
VF_iter2 = int(np.load('./Output/input_variables.npy')[3])
weight_f = np.load('./Output/input_variables.npy')[4]
weight_column_sum = np.load('./Output/input_variables.npy')[5]
VF_relax = eval(np.load('./Output/input_variables.npy')[6])
VF_stable = eval(np.load('./Output/input_variables.npy')[7])
VF_asymp = int(np.load('./Output/input_variables.npy')[8])
Passivity_Enforcement_Enable = eval(np.load('./Output/input_variables.npy')[9])
SMP_iter_upper_limit = int(np.load('./Output/input_variables.npy')[10])



print('\n********************************* Starting part one (Full synthesis) *********************************\n')

input_full_synthesis = './Input/'+ Input_file_name
TL = rf.Network(input_full_synthesis)
S_par_original = TL.s.copy()
Y_par_original = TL.y.copy()
freq = TL.frequency.f
s = 2*np.pi*1j*freq
z0 = int(TL.z0[0,0].real)
number_of_ports = S_par_original.shape[1]


# Check/Enforce reciprocity by simmetrizing S-parameters (input file)
if S_par_original.T.all() == S_par_original.all():
    S_par_reciprocal = S_par_original.copy()
    Y_par_reciprocal = Y_par_original.copy()
else:
    S_par_reciprocal = S_par_original.copy()
    # Enforce reciprocity 
    for f in range(len(freq)):
        for i in range(number_of_ports):
            for j in range(number_of_ports):
                if j<i:
                    S_par_reciprocal[f,i,j] = S_par_reciprocal[f,j,i]
    Y_par_reciprocal = admittance_from_scattering(S_par_reciprocal, z0)
    
bigS = S_par_reciprocal.copy()
bigY = Y_par_reciprocal.copy()

#---------------------------------------------------Prepare the reciprocal network (i.e., construct and stack y') for vectfit algorithm --------------------------------------- 
p = number_of_ports
Num_subckt = int((p*(p+1))/2)                                           # Num_subckt is equal to the sum of first p natural numbers (Num_subckt= p+ p-1 + ... + 1)

# Calculate Y', reshape and stack 
bigY_prime = convert_y_to_y_prime(bigY)
f_stack = np.zeros(shape=(Num_subckt, len(freq)), dtype='complex')
nt = 0
for i in range(1,p+1):
    for j in range(i,p+1):
            f_stack[nt,:] = bigY_prime[:, i-1,j-1]
            nt = nt+1

#---------------------------------------------------------------------Perform vectfit algorithm----------------------------------------------------------------------- 
f = f_stack.copy()
N = Num_poles                               # order of approximation
Ns = len(freq)
Nc = Num_subckt

Niter1 = VF_iter1                          # Fitting column sum: n.o. iterations
Niter2 = VF_iter2                          # Fitting column: n.o. iterations

#weight_f = 'sqrt_abs'                     # Choose between 'unweighted', 'abs', and 'sqrt_abs' weight types for f                            Note: 'abs' is stronger inverse wieght in comparison with 'sqrt_abs'
#weight_column_sum = 'norm'                # Choose between 'unweighted', 'norm', and 'sqrt_norm' weight types for g (column sum)             Note: 'norm' is stronger inverse wieght in comparison with 'sqrt_norm'

print('Vector Fitting: \n')

# Fitting options
opts = {}
opts['relax']     = VF_relax    # Use vector fitting with relaxed non-triviality constraint
opts['stable']    = VF_stable   # Enforce stable poles
opts['asymp']     = VF_asymp    # Fitting includes D  (asymp=0  fits with None, asymp=1  fits with D, asymp=2  fits with D and E)
opts['spy1']      = False
opts['spy2']      = False
opts['logx']      = False
opts['logy']      = False
opts['errplot']   = False
opts['phaseplot'] = False
opts['skip_pole'] = False
opts['skip_res']  = True
opts['cmplx_ss']  = True        # Will generate state space model with diagonal A
opts['legend']    = False


# Forming weight matrix for f:
if weight_f == 'unweighted':
    weight = np.ones((1,Ns))
elif  weight_f == 'abs':
    weight = 1 / np.abs(f)
elif  weight_f == 'sqrt_abs':
    weight = 1 / np.sqrt(np.abs(f))  
    

# Forming (weighted) column sum:
if weight_column_sum == 'unweighted':
    g = np.sum(f, axis=0)[nax,:]
elif  weight_column_sum == 'norm':
    g = np.sum(f / np.linalg.norm(f, axis=1)[:,nax], axis=0)[nax,:]
elif  weight_column_sum == 'sqrt_norm':
    g = np.sum(f / np.sqrt(np.linalg.norm(f, axis=1)[:,nax]), axis=0)[nax,:]
weight_g = 1 / np.abs(g)


# Complex starting poles:
bet = np.linspace(s[0].imag, s[-1].imag, int(N/2))
alf = -bet * 1.0e-2
poles = np.concatenate(((alf-1j*bet)[:,nax],(alf+1j*bet)[:,nax]), axis=1).flatten()


print('****Calculating improved initial poles by fitting column sum ...')
for it in range(Niter1):
    print('   Iter '+str(it))
    if it==Niter1-1:
        opts['skip_res'] = False
        opts['spy2']     = False
    SER,poles,rmserr,fit = vectfit3(g,s,poles,weight_g,opts)


print('****Fitting column (i.e., f) ...')
opts['skip_res'] = True
opts['spy2']     = False
for it in range(Niter2):
    print('   Iter '+str(it))
    if it==Niter2-1:
        opts['skip_res'] = False
        opts['spy2']     = False
    SER, poles, rmserr, fit = vectfit3(f,s,poles,weight,opts)

# Transforming model of lower matrix triangle into state-space model of full matrix
SER = tri2full(SER)

if opts['asymp'] == 0:
    SER['D'] = np.zeros(shape=(p,p))
    SER['E'] = np.zeros(shape=(p,p))
elif opts['asymp'] == 1:
    SER['E'] = np.zeros(shape=(p,p))
    

# Generating pole-residue model
R,a,D,E = ss2pr(SER, tri=True)
D = D.real
E = E.real
SER['R'] = R
SER['poles'] = a

print('\nEnd Vector Fitting')
print('\nRMSE = ', rmserr)
print('\n')

#------------------------------------------------Generate fitted network and compare the results with original network ----------------------------------------- 
# Genrate fitted network 
Y_fit = np.zeros(shape=(Ns,p,p), dtype='complex')
for fre in range(Ns):
    u = create_upper_matrix(fit[:,fre], p)
    l = np.triu(u,1).T
    Y_fit[fre,:,:] = u + l
    # convert from y' to y
    for i in range(p):
        Y_fit[fre,i,i] = np.sum(Y_fit[fre,i,:])
                
    for i in range(p):
        for j in range(p):
            if i != j:
                Y_fit[fre,i,j] = -Y_fit[fre,i,j]


S_par_fitted = scattering_from_admittance(Y_fit, z0)
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
    ax.semilogy(freq/1e12, np.absolute(actual_Z), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_Z), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Z)*180/np.pi, color='blue', label='Original', linestyle='solid')
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
    save_fig_VF(f"Z%d" %mn)
    
    
    # Compare Y-parameters
    actual_Y =  Y_par_reciprocal[:,m-1,n-1]
    fitted_Y = Y_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.semilogy(freq/1e12, np.absolute(actual_Y), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_Y), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Y)*180/np.pi, color='blue', label='Original', linestyle='solid')
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
    save_fig_VF(f"Y%d" %mn)
    
    
    # Compare S-parameters
    actual_S =  S_par_reciprocal[:,m-1,n-1]
    fitted_S = S_par_fitted[:,m-1,n-1]
    
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.semilogy(freq/1e12, np.absolute(actual_S), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_S), color='red', label='Fitted', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_S)*180/np.pi, color='blue', label='Original', linestyle='solid')
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
    save_fig_VF(f"S%d" %mn)
        

# Comparison plots 
for m in range(1,p+1):
    for n in range(m,p+1):
        compare_fitted_vs_actual(m,n)
        plt.close('all')

print('\n')
    
#----------------------------------------------------------------Passivity assessment/enforcement-------------------------------------------------------- 
def compare_passivefitted_vs_actual(m, n):
    mn = m*10+n
    # Compare Z-parameters
    actual_Z =  Z_par_reciprocal[:,m-1,n-1]
    fitted_Z = Z_passive_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.semilogy(freq/1e12, np.absolute(actual_Z), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_Z), color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Z)*180/np.pi, color='blue', label='Original', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_Z)*180/np.pi, color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"Z%d" %mn)
    save_fig_SMP(f"Z%d" %mn)
    
    
    # Compare Y-parameters
    actual_Y =  Y_par_reciprocal[:,m-1,n-1]
    fitted_Y = Y_passive_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.semilogy(freq/1e12, np.absolute(actual_Y), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_Y), color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_Y)*180/np.pi, color='blue', label='Original', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_Y)*180/np.pi, color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"Y%d" %mn)
    save_fig_SMP(f"Y%d" %mn)
    
    
    # Compare S-parameters
    actual_S =  S_par_reciprocal[:,m-1,n-1]
    fitted_S = S_passive_fit[:,m-1,n-1]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.semilogy(freq/1e12, np.absolute(actual_S), color='blue', label='Original', linestyle='solid')
    ax.semilogy(freq/1e12, np.absolute(fitted_S), color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5, 5)))
    ax2.set_title('Phase (deg)')
    ax2.plot(freq/1e12, np.angle(actual_S)*180/np.pi, color='blue', label='Original', linestyle='solid')
    ax2.plot(freq/1e12, np.angle(fitted_S)*180/np.pi, color='red', label='Fitted with Passivity Enforcement', linestyle=(0,(5,5)))
    ax2.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlim([freq[0]/1e12,freq[-1]/1e12])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"S%d" %mn)
    save_fig_SMP(f"S%d" %mn)
    
    
if Passivity_Enforcement_Enable == True:
    # Plot eigenvalues before passivity enforcement    
    prev_SER = SER.copy()
    ylim1 = None
    ylim2 = None
    plot(plot_name='Original', s_pass=s.T, ylim=np.array((ylim1, ylim2)), labels=['Original'], SER1=prev_SER)
    
    # Passivity Enforcement
    smp = SMP(plot=False)
    SMP_SER = smp.SMP_driver(SER, Niter=SMP_iter_upper_limit, s_pass=s.T)
    
    plot(plot_name='Orig Vs. SMP', s_pass=s.T, ylim=np.array((ylim1, ylim2)), labels=['Original', 'SMP Perturbed'], SER1=prev_SER, SER2=SMP_SER)
    
    A = SMP_SER['A']
    B = SMP_SER['B']
    C = SMP_SER['C']
    D = SMP_SER['D'].real
    E = SMP_SER['E'].real
    
    poles = SMP_SER['poles']
    residues = SMP_SER['C']
    N = poles.shape[0]
    N_port = int(residues.shape[1] / N)
    poles = poles.reshape((1, -1))
    residues = residues.reshape((N_port  ** 2, N))
    
    R = residues.reshape((N_port, N_port, N))
    residues_stack = np.zeros(shape=(Num_subckt, N), dtype='complex')
    D_stack = np.zeros(Num_subckt)
    E_stack = np.zeros(Num_subckt)
    nt = 0
    for i in range(0,p):
        for j in range(i,p):
                residues_stack[nt,:] = R[i,j,:]
                D_stack[nt] = D[i,j]
                E_stack[nt] = E[i,j]
                nt = nt+1
    
    
    Y_passive_fit = np.zeros(shape=(Ns,p,p), dtype='complex')  
    nt =0
    for s0 in s:
        A_inv = np.linalg.inv(s0*np.eye(len(A))-A)
        Y_passive_fit[nt,:,:] = C @ A_inv @ B + D
        nt = nt+1
    
    # convert from y' to y
    for fre in range(Ns):
        for i in range(p):
            Y_passive_fit[fre,i,i] = np.sum(Y_passive_fit[fre,i,:])
                    
        for i in range(p):
            for j in range(p):
                if i != j:
                    Y_passive_fit[fre,i,j] = -Y_passive_fit[fre,i,j]
    
    
    S_passive_fit = scattering_from_admittance(Y_passive_fit, z0)
    Z_passive_fit = np.linalg.inv(Y_passive_fit)
    

    # Comparison plots 
    for m in range(1,p+1):
        for n in range(m,p+1):
            compare_passivefitted_vs_actual(m,n)
            plt.close('all')
            
else:
    print('\n Passivity enforcement is disabled \n')
    
    A = SER['A']
    B = SER['B']
    C = SER['C']
    D = SER['D'].real
    E = SER['E'].real
    
    poles = SER['poles']
    residues = SER['C']
    N = poles.shape[0]
    N_port = int(residues.shape[1] / N)
    poles = poles.reshape((1, -1))
    residues = residues.reshape((N_port  ** 2, N))
    
    R = residues.reshape((N_port, N_port, N))
    residues_stack = np.zeros(shape=(Num_subckt, N), dtype='complex')
    D_stack = np.zeros(Num_subckt)
    E_stack = np.zeros(Num_subckt)
    nt = 0
    for i in range(0,p):
        for j in range(i,p):
                residues_stack[nt,:] = R[i,j,:]
                D_stack[nt] = D[i,j]
                E_stack[nt] = E[i,j]
                nt = nt+1    
    
#-----------------------------------------------------------------------Generate netlist------------------------------------------------------------------ 
print('\nGenerating Netlist:')
common_poles = np.zeros(shape=(Num_subckt, N), dtype='complex')
for i in range(Num_subckt):
    common_poles[i,:] = poles
    
create_netlist_file(common_poles, residues_stack, number_of_ports, D)


print('\n********************************* Part one (Full synthesis) is done *********************************\n')
