
""" create_netlist.py     => This script creates a netlist for a multiport network using the input poles, residues, and D matrix of fitted network generated from Run_full_synthesis.py

Author: Rasul Choupanzadeh 
Date: 08/11/2022

# Acknowledgement 1:
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [5].

Acknowledgement 2: This script is a duplication of "create_netlist" program with some 
                 modifications and corrections added by author to the original script 
                (https://github.com/JenniferEHoule/Circuit_Synthesis). All credits go to
                "Jennifer Houle" for her Python implemented program, and the following papers.

Modifications: 
  1. Used numpy.count_nonzero instead of scipy.count_nonzero in "reduce_poles" function in Jennifer's program.  Note: scipy.count_nonzero is deprecated and will be removed in SciPy 2.0.0.
  2. Corrected the Rbbr formula by adding a " - " sign in Rbbr formula in "create_imag_netlist_branch" function. 
  3. Since the shapes of input matrices differ from the input matrices of original program of "create_netlist", we did some modifications in "create_netlist_file" function.
  4. We added the option input matrx D for more accurate results. 


[1]  A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization through
     equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119

[2]  V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from Frequency Response of Multi-port
     Networks", EMC EUR, pp. 79-84, 2016. [Online]. Available: ://WOS:000392194100012.

[3]  R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent 
     circuits of improper rational transfer functions with real poles and residues. IEEE Access, vol.8, pp. 125149-125162, 2020.

[4]  Houle, Jennifer, GitHub. May 10, 2020. Accessed on: August 11, 2022, [Online]. https://github.com/JenniferEHoule/Circuit_Synthesis

[5] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

"""


## Input: poles, residues, number_of_ports, matrix D           Output: full_netlist.sp file


import numpy as np

from pathlib import Path


def reduce_poles(poles, residues):
    number_of_imaginary_poles = np.count_nonzero(poles[0, :].imag)
    for pole in range(0, poles.shape[1] - number_of_imaginary_poles // 2):
        if poles[0, pole].imag:
            poles = np.delete(poles, pole + 1, 1)
            residues = np.delete(residues, pole + 1, 1)
    return poles, residues


def create_real_netlist_branch(pole, residue, port_a, port_b, branch_number, netlist_file):
    # Equations are based on the equations in [1,3]
    print(f"* Branch {branch_number}", file=netlist_file)
    print(f"R1br{branch_number} {port_a} "
          f"net{branch_number} {(-pole.real / residue.real)}", file=netlist_file)
    print(f"L1br{branch_number} net{branch_number} "
           f"{port_b} {1 / residue.real}\n", file=netlist_file)


def create_imag_netlist_branch(pole, residue, port_a, port_b, branch_number, netlist_file):
    # Equations are based on the equations in [1]
    print(f"* Branch {branch_number}", file=netlist_file)
    print(f"Rabr{branch_number} {port_a} "
          f"netRa{branch_number} "
          f"{(residue.imag * pole.imag - residue.real * pole.real) / (2 * (residue.real) ** 2)}", file=netlist_file)
    print(f"Lbr{branch_number} netRa{branch_number} "
          f"netL{branch_number} "
          f"{1 / (2 * residue.real)}", file=netlist_file)
    print(f"Rbbr{branch_number} netL{branch_number} "
          f"{port_b} "
          f"{-((pole.imag) ** 2 * ((residue.imag) ** 2 + (residue.real) ** 2)) / (2 * (residue.real) ** 2 * (residue.imag * pole.imag + residue.real * pole.real))}", file=netlist_file)
    print(f"Cbr{branch_number} netL{branch_number} "
          f"{port_b} "
          f"{(2 * (residue.real) ** 3) / ((pole.imag) ** 2 * ((residue.imag) ** 2 + (residue.real) ** 2))}\n", file=netlist_file)



def create_netlist_branch(pole, residue, port_a, port_b, branch_number, netlist_file):
    if pole.imag:
        create_imag_netlist_branch(pole, residue, port_a, port_b, branch_number, netlist_file)
    else:
        create_real_netlist_branch(pole, residue, port_a, port_b, branch_number, netlist_file)



def create_netlist_file(poles, residues, number_of_ports, D, out_file_path='./Output/full_netlist.sp'):
    outfile = Path(out_file_path)
    with outfile.open('w') as netlist_file:
        print(f"* netlist generated with vector fitting poles and residues\n", file=netlist_file)
    
        # Main circuit declaration (configuration based on Fig. 2 in [2])
        main_declaration = ".subckt total_network "
        for port_1 in range(1, number_of_ports + 1):
            main_declaration += f"node_{port_1} "
        print(main_declaration, file=netlist_file)
    
        for port_1 in range(1, number_of_ports + 1):
            for port_2 in range(1, number_of_ports + 1):
                if port_1 == port_2:
                    port_2_name = '0'
                else:
                    port_2_name = f'node_{port_2}'
                if port_1 <= port_2:
                    print(f"X_{port_1}{port_2} node_{port_1} {port_2_name} yp{port_1}{port_2}", file=netlist_file)
        print(f".ends\n\n", file=netlist_file)     # End of the main circuit declaration
        
        # Subcircuits         
        row = 0   
        for port_1 in range(1, number_of_ports + 1):
            for port_2 in range(1, number_of_ports + 1):
                if port_1 == port_2:
                    port_2_name = '0'
                else:
                    port_2_name = f'node_{port_2}'
                if port_1 <= port_2:
                    print(f"* Y'{port_1}{port_2}", file=netlist_file)
                    print(f".subckt yp{port_1}{port_2} node_{port_1} {port_2_name}", file=netlist_file)
                    
                    poles_row = poles[row,:]
                    poles_row = poles_row.reshape(1,len(poles_row))
                    residues_row = residues[row,:]
                    residues_row = residues_row.reshape(1,len(residues_row))
                    poles_reduced, residues_reduced = reduce_poles(poles_row, residues_row)
                    for pole_num in range(poles_reduced.shape[1]):
                        create_netlist_branch(poles_reduced[0, pole_num], residues_reduced[0, pole_num], f'node_{port_1}', port_2_name, pole_num, netlist_file)
                        
                    ## Added by Rasul Choupanzadeh   
                    # Add branch for D elements 
                    if D[port_1-1,port_2-1].real != 0: 
                        branch_number_D = int(poles.shape[1]/2)
                        port_a = f'node_{port_1}'
                        port_b = port_2_name
                        print(f"* Branch {branch_number_D}", file=netlist_file)
                        print(f"Rd {port_a} "
                              f"{port_b} "
                              f"{1 / (D[port_1-1,port_2-1].real)}\n", file=netlist_file)    
                        
                    print(f".ends\n\n", file=netlist_file)
                    row = row + 1
        print(f".end", file=netlist_file)

    print(f"Wrote netlist to: {outfile.absolute()}")

