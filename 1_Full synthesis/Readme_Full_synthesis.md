# Full-order synthesis of SROPEE 
### Author: Rasul Choupanzadeh
### Date: 06/28/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program gets the optical network S-parameter results (e.g., OIDT S-parameter results in [2]) as the input, implements the vector fitting algorithm described in [3-6], and generates a full-order netlist as the output [7-10]. 

- The vector fitting part:
    * Constructs a network for the given S-parameters using scikit-rf [11]
    * Enforces the reciprocity of created network by symmetrizing S-parameeter matrix using upper triangular values of matrix
    * Calculates the Y-parameters of network
    * Fitts the Y-parameters of network using vector fitting algorithm
    * Genrates the fitted network for comparison of fitted results with the results of original network

- The netlist creation part generates a netlist (full_netlist.sp) for equivalent full-order network of given S-parameter

- The main program implements the previous two parts (i.e., vector fitting and netlist creation), and then loads the simulated results of full-order netlist from Keysight ADS software for verification of generated netlist

- Finally, it provides:
    * RMSE between fitted results and original results
    * Verification of fitted results by comparing network parameters of fitted vs. actual (original) network
    * Verification of generated full-order netlist by comparing the simulation results of netlist in Keysight ADS software vs. actual results of given S-parameter file


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Run_main_full_synthesis.py:** It is the main program file, creates the network of a given S-parameter file, calculates the Y-parameters, fitts the responses by vector fitting algorithm and provides RMSE, creates a fitted network, generates a full-order netlist file representing the equivalent network of the given network response.
    * Inputs: input file (e.g., te_rough_example.s4p), Num_pole_pairs, m,n
    * Outputs: poles, residues, full_netlist.sp file, RMSE of fitted and actual network
- **te_rough_example.s4p:** A demo input for the main program 
- **vectfit.py:** Performs the vector fitting algorithm
    * Inputs: f, s, and number of pole pairs
    * Outputs: poles and residues matrices
- **create_netlist.py:** Generates a netlist file representing the equivalent network of the given network response
    * Inputs: poles, residues, and number of number of ports
    * Outputs: full_netlist.sp file
- **full_netlist.sp:** Generated netlist file representing the equivalent circuit of the given S-parameter file
- **S_ADS_full.s4p:** Simulation results of the netlist file simulated in Keysight ADS software (note: this file is located in the folder: .\SROPEE\Output)
- **Load_ADS_full.py:** Loads the simulation results of the netlist file simulated in Keysight ADS software 
    * Inputs: S_ADS_full.s4p file
    * Outputs: S_ADS matrix



# Run instructions
To implement the full synthesis algorithm on te_rough_example.s4p input execute the Run_main_full_synthesis.py. To run the algorithm for other input file:
1. Replace your input file with existing te_rough_example.s4p file in the Full synthesis folder.
2. Replace te_rough_example.s4p with your file name in Run_main_full_synthesis.py.
2. Specify the number of pole pairs for fitting procedure.
3. You may change the m and n for comparison of different network parameters (Smn, Ymn, Zmn).


## Inputs:
- **te_rough_example.s4p:** A demo input of optical interconnects frequency response. This file may be replaced by another touchstone file of interest.
- **Num_pole_pairs:** A scalar value for number of pairs of fitting poles
- **m and n:** Scalar values for comparison of different network parameters (Smn, Ymn, Zmn)

    
## Outputs:
- **poles:** Vector of fitted poles
    * Dimensions: Number of sub-circuits (upper triangular elements of network parameter) * Number of pairs of fitting poles 
- **residues:** Vector of fitted residues
    * Dimensions: Number of sub-circuits * Number of pairs of fitting poles 
- **RMSE:** Vector of Root Mean Square Errors between the fitted Y-parameters and original Y-parameters.
    * Dimensions: 1* Number of sub-circuits 
- **full_netlist.sp:** Generated netlist file for equivalent circuit of given S-parameter file.



# Usage:
This program was designed to provide an equivalent netlist for the frequency responses of optical interconnects, which is the initial step of SROPEE before order reductio.

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

[9] R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent circuits of improper rational transfer functions with real poles and residues. IEEE Access, vol.8, page.125149:125162, 2020.

[10] https://github.com/JenniferEHoule/Circuit_Synthesis

[11] http://scikit-rf.org

```
