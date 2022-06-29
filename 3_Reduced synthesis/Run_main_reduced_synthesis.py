
""" Run_main_reduced_synthesis.py      => This program implements the reduced circuit synthesis for the given reduced MNA matrices, and exports the equivalent reduced netlist as reduced_netlist.sp file

Author: Rasul Choupanzadeh
Date: 06/28/2022

"""


import numpy as np
import matplotlib.pyplot as plt
from Load_equiv_circuit_results import Load_equiv_circuit_results

# IMPORTANT NOTICE : This program is the continuation of passive MOR part of SROPEE. For execution of full program, you should: 
#                           1. Load the reduced data using reduced_MNA_data.npz file, which is generated in the passive MOR part of SROPEE
#                           2. Run the line "exec(open("equivalent_circuit.py").read())" to generate "equiv_netlist.sp"     
#                           2. Imoprt the "equiv_netlist.sp" to ADS and export the result of ADS as "reduced_voltages.cti"
#                           3. Run from line "V_equiv = Load_equiv_circuit_results(n,Q)" to the end, which loads the results of ADS of reduced equivalent circuit for comparison
#
# Therefore, without performing step 3, some lines will not be run or will run according to the previous generated "reduced_voltages.cti" file, which belongs to the case of input current only in port 1 (i.e., Z11, Z21, Z31, Z41)



## inputs : reduced_MNA_matrices.npz            output : reduced_netlist.sp file



# -------------------------------------------------Load the reduced Matrices and data----------------------------------------------------
data = np.load('reduced_MNA_data.npz')
Cmr = data['name1']
Gmr = data['name2']
Gammar = data['name3']
Bmr = data['name4']
Q = data['name5']
Jm = data['name6']
in_port = data['name7']
out_port = data['name8']
freq = data['name9']

# -------------------------------------------------Run Reverse MNA Algorithm----------------------------------------------------
exec(open("equivalent_circuit.py").read())

# -----------------------------------Load the simulation results of generated reduced netlist-----------------------------------
## This part Loads the pre-generated data for a specific network, therefore, you should generate your own data and replace it with "reduced_voltages.cti".
V_equiv = Load_equiv_circuit_results(n,Q)
V_ports_equiv = V_equiv[0:p,:]
Z_equiv = V_ports_equiv[out_port,:]/Jm[in_port]
plt.semilogy(freq, np.absolute(Z_equiv), label='Equivalent circuit', c="green", linestyle='dotted')
plt.legend(loc='upper left')
plt.xlim([freq[0],freq[-1]])
plt.xlabel('Frequency (THz)')
plt.ylabel('$|Z|$')
plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)

plt.show()