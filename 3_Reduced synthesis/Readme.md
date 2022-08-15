# Reduced-order synthesis part of SROPEE algorithm
### Author: Rasul Choupanzadeh
### Date: 08/15/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program generates a reduced-order equivalent netlist for the exported reduced-order MNA matrices from Passive MOR part of SROPEE.
- It loads the reduced-order MNA matrices/data from the exported file from passive MOR part of SROPEE.
- It runs the Reverse MNA algorithm to implement the circuit synthesis of reduced-order system. This program:
    * Calculates the circuit parameters of equivalent reduced-order circuit
    * Generates a reduced-order netlist (reduced_netlist.sp) representing the equivalent reduced-order circuit.
- Finally, it loads the simulation results of reduced-order netlist, simulated in Keysight ADS software, for the purpose of verification  



# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **reduced_MNA_data.npz:** Reduced-order MNA matricres and data generated from passive MOR part of SROPEE (note: this file is located in the folder: .\SROPEE\Output).
- **equivalent_circuit.py:** Calculates the parameter of equivalent reduced-order circuit, and generates the equivalent reduced-order netlist as reduced_netlist.sp file.
    * Inputs: Cmr,Gmr,Gammar,Bmr,Jm
    * Outputs: reduced_netlist.sp 
- **reduced_netlist.sp:** Generated equivalent reduced-order netlist due to the given in_port.
- **reduced_voltages_ADS.cti:** Simulation results of generated reduced-order netlist in Keysight ADS software (note: this file is located in the folder: .\SROPEE\Output).
- **Load_equiv_circuit_results.py:** Loads the simulation results of reduced-order netlist from reduced_voltages_ADS.cti. 
    * Inputs: n, Q, and "reduced_voltages_ADS.cti" file
    * Outputs: V_equiv
- **Run_reduced_synthesis.py:** This is the main program file for Reduced synthesis part of SROPEE. It loads the reduced MNA matrices and data, calculates the parameters of reduced-order circuit, generates the equivalent reduced-order netlist, and loads the simulation results of reduced-order netlist simulated in Keysight ADS.



# Run instructions
This program is a part of SROPEE algorithm, and will be executed auromatically by executing the main program file of SROPEE (run_main.py)


## Inputs:
- **in_port:** An input variable specifying the input port for generating related reduced-order netlist and calculation of reduced-order nodal voltages in Keysight ADS software (note: this variable is loaded from input_variables.npy located in the folder: .\SROPEE\Output)
- **reduced_MNA_data.npz:** Reduced-order MNA matricres and data generated from passive MOR part of SROPEE (note: this file is located in the folder: .\SROPEE\Output)
- **reduced_voltages_ADS.cti:** Reduced-order nodal voltages from simulation of reduced-order netlist in Keysight ADS software, used for visualization/verification of Z-parameter (note: this file is located in the folder: .\SROPEE\Outputnote, which is generated with the assumption of current source only in port 1, therefore, it will be used for plotting and comparing Z-parameters with inputs in port 1, (i.e, Z11, Z21, Z31, Z41). It may be replaced by another *.cti files generated with different assumptions)


    
## Outputs:
- **reduced_netlist.sp:** Generated equivalent reduced-order netlist (note: this netlist is generated according the in_port variable, which is shown as .param Ip=1.0 for that port, for the other assumptions of input ports, you need only set the related .param Ip=1.0 for that port, instead of generating a new reduced-order netlist)
- **figures:** Plots of Z-parameters of network, in magnitude and phase, according to obtained reduced-order nodal voltages form simulation of equivalent reduced-order netlist in Keysight ADS software (note: The graphs are saved as .pdf files in the folder: .\SROPEE\Output\figures\Reduced_synthesis)



# Usage:
This program was designed to provide an equivalent reduced-order netlist for a given optical interconnects. It is the final part (third part) of SROPEE.

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
