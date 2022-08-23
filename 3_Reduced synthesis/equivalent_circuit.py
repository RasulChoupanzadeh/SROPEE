
""" equivalent_circuit.py      => This script calculates parameters of equivalent circuit for reduced-order system, then, generates the reduced-order netlist as "reduced_netlist.sp"

Author: Rasul Choupanzadeh
Date: 07/03/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

"""


## inputs : Cmr,Gmr,Gammar,Bmr,Jm             output : reduced_netlist.sp file


import numpy as np

p = Bmr.shape[1]
n = len(Cmr)

# Calculating parameters of equivalent circuit
Ceq = -Cmr
Geq = -Gmr
Gammaeq = -Gammar
for i in range(0,n):
    Ceq[i,i] = -np.sum(Ceq[i,:])
    Geq[i,i] = -np.sum(Geq[i,:])
    Gammaeq[i,i] = -np.sum(Gammaeq[i,:])
    
Req = 1/Geq[:]
Leq = 1/Gammaeq[:]
Ieq = Bmr

# Converting upper triangular values of equivalent matrices into a vector to write into the equivalent netlist
R = Req[np.triu_indices(n)]
L = Leq[np.triu_indices(n)]
C = Ceq[np.triu_indices(n)]

# Writing the names and values in equivalent netlist file----------------------------------
with open('./Output/reduced_netlist.sp', 'w') as f:
    f.write('* netlist generated with reverse MNA (number of voltage nodes: n = ' + str(n) + ' )\n\n')
    f.write('.subckt equivalent_circuit\n\n')
    
    # Writing parameters I1, I2,...,Ip
    for i in range(p):
        f.write('.param Ip' + str(i+1) + '=' + str(float(Jm[i])) + '\n')
    f.write('\n')
    
    # Writing R, L, C values
    nt = 0
    for i in range(n):
        for j in range(i,n):
            n1 = 'V' + str(i+1)
            if j==i:
                n2 = '0'
            else:
                n2 = 'V' + str(j+1)
            f.write('R' + str(i+1) + '_' + str(j+1) + ' ' + n1 + ' ' + n2 + ' ' + str(float(R[nt])) + '\n')
            f.write('L' + str(i+1) + '_' + str(j+1) + ' ' + n1 + ' ' + n2 + ' ' + str(float(L[nt])) + '\n')
            f.write('C' + str(i+1) + '_' + str(j+1) + ' ' + n1 + ' ' + n2 + ' ' + str(float(C[nt])) + '\n\n')
            nt = nt + 1
            
    # Writing currnet values
    for i in range(n):
        n1 = '0'
        n2 = 'V' + str(i+1)
        for j in range(p):
            f.write('ISRC' + str(i+1) + '_p' + str(j+1) + ' ' + n1 + ' ' + n2 + ' ' + ' ac' + ' ' + str("{:0.15e}".format(Ieq[i,j])) + '*Ip' + str(j+1) + '\n')
        f.write('\n')
        
    # Closing notations of equivalent netlist
    f.write('\n' + '.ends subckt equivalent_circuit\n\n')
    f.write('.end\n\n')
        