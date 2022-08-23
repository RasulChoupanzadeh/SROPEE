
""" Run_passive_MOR.py      => This script is the main program of passive model order reduction (Passive MOR) part of SROPEE. 
Author: Rasul Choupanzadeh
Date: 08/14/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [6].

In addition to Load_netlist.py, this code uses general_MNA_builder.py and Block_SAPOR.py (and SOrth.py), which are based on the concepts from [1-2] and [3-4], respectively.

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6), pp. 504-509, 1975.
    
[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

[3] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pp. 74-79, 2004.
    
[4] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pp. 244-249, 2005.
    
[5] http://scikit-rf.org

[6] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

"""



## Input: freq (from Input_file_name), n, full_netlist.sp          Output: reduced_MNA_matrices.npz, comparison figures between full-order and reduced-order results


import os
import numpy as np
import matplotlib.pyplot as plt
import time
from general_MNA_builder import MNA_calculator
from sklearn.metrics import mean_squared_error
from math import sqrt
import skrf as rf
from skrf import Network, Frequency
      
        
# Where to save the figures
PROJECT_ROOT_DIR = "./Output/"
CHAPTER_ID = "Passive_MOR"
IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "figures", CHAPTER_ID)
os.makedirs(IMAGES_PATH, exist_ok=True)


def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)    

        
#----------------------------------------------------------------------Input parameters-----------------------------------------------------------------------------------
Input_file_name = np.load('./Output/input_variables.npy')[0]
n = int(np.load('./Output/input_variables.npy')[11])

print('\n********************************** Starting part two (Passive MOR) **********************************\n')

# creating frequency matrix 
input_passive_MOR = './Input/'+ Input_file_name
TL = Network(input_passive_MOR)   
freq_start = TL.frequency.start   
freq_stop = TL.frequency.stop
Num_freq = TL.frequency.npoints
freq = np.linspace(freq_start,freq_stop,Num_freq)

Coef = 1e12  
freq = freq/Coef
fmin = freq[0]
fmax = freq[-1]
#s0 = 2*np.pi*float((fmin+fmax)/2)*Coef  
s0 = 100*Coef

#-----------------------------------------------------------------


#load the netlist
exec(open("./2_Passive MOR/Load_netlist.py").read()) 
                                                      
#-----------------------------------------------------------------General MNA approach for MOR-------------------------------------------------------------------------
# build general MNA
[Gm,Cm,Gamma,Lv,Bm,M,N,SEt] = MNA_calculator(Ra,Rb,L,C,Rd,Num_branch, Num_port);
p = Bm.shape[1];
q = Lv.shape[1];

# Run Block SAPOR Algorithm for Model Order Reduction and obtain reduced matrices
exec(open("./2_Passive MOR/Block_SAPOR.py").read())

np.savez('./Output/reduced_MNA_data.npz', name1=Cmr, name2=Gmr, name3=Gammar, name4=Bmr, name5=Q, name6=freq)
print('Reduced MNA matrix data is saved as reduced_MNA_data.npz in output folder')


#-----------------------------------------------------------Frequency response calculation-----------------------------------------------------------------------------
## If you want to speed up the SROPEE, you can disable this whole or part of (full system part) this section (i.e., Frequency response calculation)
RMSE = []
t_full = 0;
t_reduced = 0;
for m in range(1,p+1):
    for n in range(m,p+1):
        in_port = n-1
        out_port = m-1
        Jm = np.zeros(shape=(Num_port,1)) 
        Jm[in_port] = 1 

        # Frequency response of full system
        start = time.time()
        V_orig = np.zeros(shape=(q,len(freq)), dtype = 'complex')
        nt = 0
        for f in freq:        
            s = 1j*2*np.pi*f*Coef
            V = np.transpose(Lv) @ (np.linalg.inv(s*Cm + Gm + Gamma/s) @ Bm) @ Jm
            V_orig[:,nt]=V[:,0]
            nt = nt+1
        end = time.time()              
        t_full = t_full + end-start;
        V_ports_full = V_orig[0:p,:]
        Z_MNA_full = V_ports_full[out_port,:]/Jm[in_port]
        
        
        # Frequency response of reduced system
        start = time.time()
        V_reduced = np.zeros(shape=(q,len(freq)), dtype = 'complex')   
        nt=0
        for f in freq:
            s = 1j*2*np.pi*f*Coef
            V = np.transpose(Lvr) @ np.linalg.inv(s*Cmr + Gmr + Gammar/s) @ Bmr @ Jm
            V_reduced[:,nt] = V[:,0]
            nt = nt+1
        end = time.time()
        t_reduced = t_reduced + end-start;
        V_ports_reduced = V_reduced[0:p,:]
        Z_MNA_reduced = V_ports_reduced[out_port,:]/Jm[in_port]
        
        #plot
        mn = m*10+n
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax.set_title('Magnitude')
        ax.plot(freq, np.absolute(Z_MNA_full), color='blue', label='MNA Full-order', linestyle='solid')
        ax.plot(freq, np.absolute(Z_MNA_reduced), color='red', label='MNA Reduced-order', linestyle=(0,(5, 5)))
        ax2.set_title('Phase (deg)')
        ax2.plot(freq, np.angle(Z_MNA_full)*180/np.pi, color='blue', label='MNA Full-order', linestyle='solid')
        ax2.plot(freq, np.angle(Z_MNA_reduced)*180/np.pi, color='red', label='MNA Reduced-order', linestyle=(0,(5,5)))
        ax2.set_xlim([fmin,fmax])
        ax.set_xlim([fmin,fmax])
        ax.set_xlabel('Frequency (THz)')
        ax2.set_xlabel('Frequency (THz)')
        ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
        ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
        ax.legend(loc='upper right')
        ax2.legend(loc='upper right')
        fig.suptitle(f"Z%d" %mn)
        save_fig(f"Z%d" %mn)
        
        # RMSE
        square = np.square(np.absolute(Z_MNA_reduced-Z_MNA_full))
        MSE=square.mean()
        rmse = np.sqrt(MSE)
        RMSE.append(rmse)        

plt.close('all')

# Info
print("The execution time of full-order MNA system is", t_full, 'seconds')
print("The execution time of reduced-order MNA system is", t_reduced, 'seconds')
print("Maximum RMSE of reduced and full order MNA = ", np.max(RMSE))



#---------------------------------------------------------------Statement---------------------------------------------------------------------------------
print('\n********************************** Part two (Passive MOR) is done ***********************************\n')
