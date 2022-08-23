
""" Run_main_reduced_synthesis.py      => This program implements the reduced circuit synthesis for the reduced MNA matrices exported from Passive MOR part of SROPEE, and exports the equivalent reduced netlist as reduced_netlist.sp file

Author: Rasul Choupanzadeh
Date: 07/03/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

"""


## inputs : in_port, reduced_MNA_matrices.npz            output : reduced_netlist.sp file


import os
import numpy as np
import matplotlib.pyplot as plt
from Load_equiv_circuit_results import Load_equiv_circuit_results

# IMPORTANT NOTICE : This program is the continuation of passive MOR part of SROPEE. For execution of full program, you should: 
#                           1. Run to the line "exec(open("equivalent_circuit.py").read())" to generate "equiv_netlist.sp" according your port of input source (n)   
#                           2. Imoprt the "equiv_netlist.sp" to ADS and export the result of ADS as "reduced_voltages.cti"
#                           3. Run from line "V_equiv = Load_equiv_circuit_results(n,Q)" to the end, which loads the results of ADS of reduced equivalent circuit for comparison
#
# Therefore, without performing step 3, some lines will not be run or will run according to the previous generated "reduced_voltages.cti" file, which belongs to the case of input current only in port 1 (i.e., Z11, Z21, Z31, Z41)


# Where to save the figures
PROJECT_ROOT_DIR = "./Output/"
CHAPTER_ID = "Reduced_synthesis"
IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "figures", CHAPTER_ID)
os.makedirs(IMAGES_PATH, exist_ok=True)


def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)    



# -------------------------------------------------Load the reduced Matrices and data-------------------------------------------
in_port = int(np.load('./Output/input_variables.npy')[12])               
 

print('\n****************************** Starting part three (Reduced synthesis) ******************************\n')

data = np.load('./Output/reduced_MNA_data.npz')
Cmr = data['name1']
Gmr = data['name2']
Gammar = data['name3']
Bmr = data['name4']
Q = data['name5']
freq = data['name6']

fmin = freq[0]
fmax = freq[-1]
Num_port = Bmr.shape[1]
Jm = np.zeros(shape=(Num_port,1)) 
Jm[in_port-1] = 1 

n = Cmr.shape[0]

# -------------------------------------------------Run Reverse MNA Algorithm----------------------------------------------------
exec(open("./3_Reduced synthesis/equivalent_circuit.py").read())
print('Equivalent reduced netlist is saved as reduced_netlist.sp in output folder')




# -----------------------------------Load the simulation results of generated reduced netlist-----------------------------------
## If you want to speed up the SROPEE, you can disable this section (i.e., Load the simulation results of generated reduced netlist)
## This part Loads the pre-generated data for a specific network with input source at port 1, therefore, you should generate your own data and replace it with "reduced_voltages.cti".
V_equiv = Load_equiv_circuit_results(n,Q)
V_ports_equiv = V_equiv[0:p,:]

for m in range(1,Num_port+1):
    out_port = m-1
    Z_equiv_reduced = V_ports_equiv[out_port,:]/Jm[in_port-1]
    
    mn = m*10+in_port
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax.set_title('Magnitude')
    ax.plot(freq, np.absolute(Z_equiv_reduced), color='blue', label='Equivalent reduced circuit', linestyle='solid')
    ax2.set_title('Phase (deg)')
    ax2.plot(freq, np.angle(Z_equiv_reduced)*180/np.pi, color='blue', label='Equivalent reduced circuit', linestyle='solid')
    ax2.set_xlim([freq[0],freq[-1]])
    ax.set_xlim([freq[0],freq[-1]])
    ax.set_xlabel('Frequency (THz)')
    ax2.set_xlabel('Frequency (THz)')
    ax.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax2.grid(color = 'green', linestyle = '--', linewidth = 0.5)
    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')
    fig.suptitle(f"Z%d" %mn)
    save_fig(f"Z%d" %mn)    


#---------------------------------------------------------------Statement---------------------------------------------------------------------------------
print('\n****************************** Part three (Reduced synthesis) is done *******************************\n')

