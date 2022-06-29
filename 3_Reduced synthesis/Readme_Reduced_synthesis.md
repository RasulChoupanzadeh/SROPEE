# Reduced-order synthesis of SROPEE 
### Author: Rasul Choupanzadeh
### Date: 06/28/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program generates a reduced equivalent netlist, representing reduced equivalent circuit, for the given reduced MNA matrices (reduced system).

- It loads the reduced data from reduced_MNA_data.npz file, which is generated in the passive MOR part of SROPEE

- It runs the Reverse MNA algorithm. The Reverse MNA part part:
    * Calculates the circuit parameters of equivalent reduced circuit
    * Generates a reduced netlist (reduced_netlist.sp) representing the equivalent reduced circuit.

- Finally, it loads the simulation results of reduced netlist, simulated in Keysight ADS software, for the purpose of verification  


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **reduced_MNA_data.npz:** Reduced MNA matricres and data generated from passive MOR part of SROPEE. This file includes Cmr,Gmr,Gammar,Bmr,Q,Jm,in_port,out_port, and freq.
- **equivalent_circuit.py:** Calculates the parameter of equivalent reduced circuit, and generates the equivalent reduced netlist as reduced_netlist.sp file.
    * Inputs: Cmr,Gmr,Gammar,Bmr,Jm
    * Outputs: reduced_netlist.sp 
- **reduced_netlist.sp:** Generated equivalent reduced netlist
- **reduced_voltages.cti:** Simulation results of reduced netlist in ADS software.
- **Load_equiv_circuit_results.py:** Loads the simulation results of reduced netlist. 
    * Inputs: n, Q, and "reduced_voltages.cti" file
    * Outputs: V_equiv
- **Run_main_reduced_synthesis.py:** It is the main program file, which loads the reduced MNA matrices and data, calculates the parameters of reduced circuit, generates the equivalent reduced netlist, and loads the simulation results of reduced netlist, by calling and/or running all the previously mentioned program files.



# Run instructions
To implement the reduced synthesis algorithm on the reduced system:
1. Copy your MNA matrices and data file (e.g., reduced_MNA_data.npz) in Reduced synthesis folder.
2. Execute Run_main_reduced_synthesis.py with a breakpoint at the end of Reverse MNA part of Run_main_reduced_synthesis.py.
2. Load the reduced_netlist.sp in ADS software, simulate, and extract the results as reduced_voltages.cti according to the instructions.pdf of SROPEE.
3. Execute the remainder of Run_main_reduced_synthesis.py.


## Inputs:
- **reduced_MNA_data.npz:** Reduced MNA matricres and data generated from passive MOR part of SROPEE.
- **reduced_voltages.cti:** Simulation results of reduced netlist in ADS software.

    
## Outputs:
- **reduced_netlist.sp:** Generated equivalent reduced netlist 
- **Figure:** Plots of simulated Z-parameter of equivalent reduced netlist in ADS software.



# Usage:
This program was designed to provide an equivalent reduced netlist for the reduced order system of optical interconnects, which is the final step (third step) of SROPEE.

# Software Version Information:
**Python 3.8.13**

Libraries used in Python:
   * pip		21.2.2
   * numpy		1.21.5
   * matplotlib	        3.5.1

# References:
```
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects," National Science Foundation, Award #1816542. Jun. 22, 2018.

```
