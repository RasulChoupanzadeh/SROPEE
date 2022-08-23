# Full-order synthesis part of SROPEE algorithm
### Author: Rasul Choupanzadeh
### Date: 08/23/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program gets the optical network S-parameter results (e.g., OIDT S-parameter results in [2]) as the input, implements the vector fitting algorithm described in [3-6], assesses/enforces the passivity by implementing algorithm described in [7,8], and generates a full-order netlist as the output [9-12]. 

- The vector fitting part:
    * Constructs a network for the given S-parameters using scikit-rf [13]
    * Enforces the reciprocity of created network by symmetrizing S-parameeter matrix using upper triangular values of matrix
    * Calculates the Y-parameters of network
    * Fits the Y-parameters of network using vector fitting algorithm [6]
    * Genrates a fitted network for comparison of fitted results with the results of original (full-order) network

- The Singular Matrix Passivity (SMP) enforcement part:
    * Completes the state space model (SER)
    * Checks the passivity of fitted model by obtaining eigenvalues
    * Enforces the passivity of non-passive fitted models [7,8]
- The netlist creation part generates a netlist (full_netlist.sp) representing equivalent circuit for full-order network of given S-parameter

- The Run_full_synthesis.py is the main program of full synthesis part of SROPEE, which implements previously mentioned vector fitting and netlist creation programs.

- Finally, it:
    * Provides RMSE between passive fitted results vs. original results
    * Provides verification of fitted results by comparing the network parameters of fitted vs. actual (original full-order) network
    * Saves the outputs in the output folder of SROPEE


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Run_full_synthesis.py:** It is the main program file for full synthesis part of SROPEE
    * Inputs: Input_file_name and Num_poles from input_variables.npy 
    * Outputs: Poles, residues, full_netlist.sp file, RMSE of fitted and actual network, comparison figures of actual vs. fitted network
- **THz_S_param.s4p:** A demo input for the main program (note: this file is located in the folder: .\SROPEE\Input)
- **input_variables.npy:** Input variables of Run_full_synthesis.py (note: this file is located in the folder: .\SROPEE\Output)
- **vectfit3.py:** This function is imported from [6], and performs the vector fitting algorithm
    * Inputs: f, s, and number of poles
    * Outputs: poles and residues matrices
- **SMP.py:** Singularity matrix perturbation implementation of [7] imported from [8]
- **pr2ss.py:** This function imported from [12], and completes the state space model for passivity assessment/enforcement.
- **intercheig.py, rot.py:** This function imported from [12], and adjusts the eigenvalues / eigenvectors. 
- **eig_plot.py:** This function is imported from [8], and provides eigenvalue plots
- **create_netlist.py:** Generates a netlist file representing the equivalent network of the given network responses
    * Inputs: poles, residues, D matrix, and number of number of ports
    * Outputs: full_netlist.sp file 
- **full_netlist.sp:** Generated netlist file representing the equivalent circuit of the given S-parameter file (note: this file is saved in the folder: .\SROPEE\Output)



# Run instructions
This program is a part of SROPEE algorithm, and will be executed auromatically by executing the main program file of SROPEE (run_main.py)


## Inputs:
- **THz_S_param.s4p:** A demo input of optical interconnects frequency response. This file is located in the input folder of SROPEE, and may be replaced by another touchstone file of interest.
- **Num_poles:** Number of fitting poles (approximation order in vector fitting)
- **VF_iter1:** Iteration number for calculating improved initial poles by fitting column sum in vector fitting algorithm
- **VF_iter2:** Iteration number for calculating poles and residues in vector fitting algorithm
- **weight_f:** Weighting type for vector fitting of function (f)
- **weight_column_sum:** Weighting type for fitting of column sum
- **VF_relax:** Enables/disables relaxed vector fitting
- **VF_stable:** Enables/disables enforcing stability of fitted poles in vector fitting algorithm
- **VF_asymp:** Fitting model options (fit with None, fit with D, fit with D and E)
- **Passivity_Enforcement_Enable:** Enables/disables passivity assessment/enforcement
- **SMP_iter_upper_limit:** Passivity enforcement iteration upper limit
- **input_variables.npy:** This file contains above inputs, and is used here to load the mentioned inputs (note: this file is located in the folder: .\SROPEE\Output)

    
## Outputs:
- **SER:** State space model in a dictionary:
    * A: Final poles in a diagonal matrix
	  * Dimension: Num_poles * Num_poles
    * B: Vector of 1's
	  * Dimension: Num_poles * 1
    * C: Residues in a matrix (each row corresponds to a single Y' frequency response)
	  * Dimension: number of subcircuits (Num_sbck) * Num_poles
    * D: d parameters in fitting model
	  * Dimension: Num_sbck * Num_poles
    * E: e parameters in fitting model
	  * Dimension: Num_sbck * Num_poles
- **poles:** Vector of fitted poles
    * Dimensions: 1* Num_poles
- **residues:** Vector of fitted residues
    * Dimensions: Number of sub-circuits * Num_poles 
- **rmserr:** RMS error from the fit generated with SER
    * Dimensions: number_of_iterations X 1
- **RMSE:** Root Mean Square Errors between passive_fitted Y-parameters and actual (original full-order) Y-parameters
- **full_netlist.sp:** Generated netlist file for equivalent circuit of given S-parameter file (note: this file is saved in the folder: .\SROPEE\Output)
- **figures:** Generated fitted vs. actual network, and passive_fitted vs. actual network graphs shown in magnitude and phase (note: The graphs are saved as .pdf files in the folder: .\SROPEE\Output\figures\Full_synthesis)


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

[2] Guiana, Brian, GitHub. May 14, 2021. Available:https://github.com/bmguiana/OIDT

[3] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

[4] B. Gustavsen, "Improving the pole relocating properties of vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

[5] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method", IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.
     
[6] Simon De Ridder, GitHub. Feb 08, 2019. Accessed on: August 10, 2022, [Online].https://github.com/SimonDeRidder/pyVF

[7] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery, vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[8] Totorica, Nathan, GitHub, May 5, 2021. Accessed on: August 10, 2022, [Online]. Available:https://github.com/ntotorica/SMP_Passivity_Enforcement

[9] A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization through equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119

[10] V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from Frequency Response of Multi-port Networks", EMC EUR, pp. 79-84, 2016. [Online]. Available: ://WOS:000392194100012.

[11] R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent circuits of improper rational transfer functions with real poles and residues. IEEE Access, vol.8, pp. 125149-125162, 2020.

[12] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021, [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis

[13] http://scikit-rf.org

```




 
