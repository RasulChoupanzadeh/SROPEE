# Full-order synthesis part of SROPEE algorithm
### Author: Rasul Choupanzadeh
### Date: 07/03/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program gets the optical network S-parameter results (e.g., OIDT S-parameter results in [2]) as the input, implements the vector fitting algorithm described in [3-6], and generates a full-order netlist as the output [7-10]. 

- The vector fitting part:
    * Constructs a network for the given S-parameters using scikit-rf [10]
    * Enforces the reciprocity of created network by symmetrizing S-parameeter matrix using upper triangular values of matrix
    * Calculates the Y-parameters of network
    * Fits the Y-parameters of network using vector fitting algorithm
    * Genrates a fitted network for comparison of fitted results with the results of original (full-order) network

- The netlist creation part generates a netlist (full_netlist.sp) representing equivalent circuit for full-order network of given S-parameter

- The Run_full_synthesis.py is the main program of full synthesis part of SROPEE, which implements previously mentioned vector fitting and netlist creation programs.

- Finally, it:
    * Provides RMSE between fitted results and original results
    * Provides verification of fitted results by comparing the network parameters of fitted vs. actual (original full-order) network
    * Saves the outputs in the output folder of SROPEE


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Run_full_synthesis.py:** It is the main program file for full synthesis part of SROPEE
    * Inputs: Input_file_name and Num_pole_pairs from input_variables.npy 
    * Outputs: Poles, residues, full_netlist.sp file, RMSE of fitted and actual network, comparison figures of actual vs. fitted network
- **te_rough_example.s4p:** A demo input for the main program (note: this file is located in the folder: .\SROPEE\Input)
- **input_variables.npy:** Input variables of Run_full_synthesis.py (note: this file is located in the folder: .\SROPEE\Output)
- **vectfit.py:** Performs the vector fitting algorithm
    * Inputs: f, s, and number of pole pairs
    * Outputs: poles and residues matrices
- **create_netlist.py:** Generates a netlist file representing the equivalent network of the given network responses
    * Inputs: poles, residues, and number of number of ports
    * Outputs: full_netlist.sp file 
- **full_netlist.sp:** Generated netlist file representing the equivalent circuit of the given S-parameter file (note: this file is saved in the folder: .\SROPEE\Output)



# Run instructions
This program is a part of SROPEE algorithm, and will be executed auromatically by executing the main program file of SROPEE (run_main.py)


## Inputs:
- **te_rough_example.s4p:** A demo input of optical interconnects frequency response. This file is located in the input folder of SROPEE, and may be replaced by another touchstone file of interest.
- **Num_pole_pairs:** A scalar value for number of pairs of fitting poles
- **input_variables.npy:** This file contains input variables, and is used here to load the input file name and Num_pole_pairs (note: this file is located in the folder: .\SROPEE\Output)

    
## Outputs:
- **poles:** Vector of fitted poles
    * Dimensions: Number of sub-circuits (upper triangular elements of network parameter) * Number of pairs of fitting poles 
- **residues:** Vector of fitted residues
    * Dimensions: Number of sub-circuits * Number of pairs of fitting poles 
- **RMSE:** Vector of Root Mean Square Errors between the fitted Y-parameters and actual (original full-order) Y-parameters
    * Dimensions: 1* Number of sub-circuits 
- **full_netlist.sp:** Generated netlist file for equivalent circuit of given S-parameter file (note: this file is saved in the folder: .\SROPEE\Output)
- **figures:** Generated fitted vs. actual network graphs shown in magnitude and phase (note: The graphs are saved as .pdf files in the folder: .\SROPEE\Output\figures\Full_synthesis)


# Usage:
This program was designed to provide an equivalent full-order netlist for the frequency responses of optical interconnects, and it is the initial step of SROPEE algorithm.

# Software Version Information:
**Python 3.8.13**

Libraries used in Python:
   * pip		21.2.2
   * scikit-learn	1.0.2
   * numpy		1.21.5
   * matplotlib	        3.5.1
   * scikit-rf          0.21.0



# References:
```
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

[2] https://github.com/bmguiana/OIDT/tree/main/CEM/GPU/FDTD_3D/s-parameter_extraction/Results

[3] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

[4] B. Gustavsen, "Improving the pole relocating properties of vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

[5] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method", IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.
     
[6] https://github.com/PhilReinhold/vectfit_python

[7] A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization through equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119

[8] V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from Frequency Response of Multi-port Networks", EMC EUR, pp. 79-84, 2016. [Online]. Available: ://WOS:000392194100012.

[9] R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent circuits of improper rational transfer functions with real poles and residues. IEEE Access, vol.8, pp. 125149-125162, 2020.

[10] https://github.com/JenniferEHoule/Circuit_Synthesis

[11] http://scikit-rf.org

```
