# S-parameter to Reduced-Order Passivity-Enforced Equivalent Circuit (SROPEE)
### Author: Rasul Choupanzadeh
### Date: 08/23/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This algorithm has three main parts in three different folders in the following order:
    1. Full synthesis
    2. Passive Model Order Reduction (Passive MOR)
    3. Reduced synthesis

- The Full synthesis part implements a full-order circuit synthesis on a given S-parameter file (input file) by performing:
    * Vector Fitting algorithm to obtain fitted poles/residues and state space model [3-6]
    * Passivity assessment/enforcement by singularity test matrix method [7,8]
    * Calculating full-order circuit parameters from poles/residues [9-11]
    * Generating full-order SPICE netlist from obtained circuit parameters [12]

- The Passvie MOR part implements Block SAPOR as the model order reduction technique on the full-order SPICE netlist file exported from Full synthesis part of SROPEE. It performs:
    * General MNA method for calculating full-order MNA matrices [14]
    * Block SAPOR algorithm to calculate reduced-order MNA matrices [16,17]

- The Reduced synthesis part implements circuit synthesis on the reduced-order MNA matrices exported from Passive MOR part. It performs:
    * Reverse MNA algorithm to obtain reduced-order circuit parameters from reduced-order MNA matrices
    * Generating reduced-order SPICE netlist from obtained reduced-order circuit parameters 

- Finally, it provides an equivalent reduced-order netlist, which may be imported to Keysight ADS software for further simulations and verification of obtained results.



# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Full synthesis folder:** Implements full-order circuit synthesis to generate equivalent full-order netlist.
    * Inputs: THz_S_param.s4p (note: this file can be replaced by any other touchstone file of interest)
    * Outputs: full_netlist.sp 
- **Passive MOR:** Implements model order reduction algorithm to calculate reduced-order MNA and projection matrices.
    * Inputs: full_netlist.sp 
    * Outputs: reduced_MNA_data.npz
- **Reduced synthesis folder:** Implements reduced-order circuit synthesis to generate equivalent reduced-order netlist.
    * Inputs: reduced_MNA_data.npz
    * Outputs: reduced_netlist.sp 
- **MNA_Methods folder:** This folder provides a comparison between Piece-by-Piece and General MNA Approaches for Model Order Reduction of Optical multiport Networks via SAPOR/Block SAPOR Algorithms. This work is done as an additional research about the proposed approaches of creating MNA matrices, and it is not a main part of the SROPEE algorithm.
- **Flow chart.pdf:** Provides the general flow chart of SROPEE algorithm.
- **Instruction.pdf:** Provides fully detailed steps for executing SROPEE algorithm. 


# Run instructions
To see the fully detailed instructions, please see the Instructions.pdf file.


## Inputs:
- **THz_S_param.s4p:** A demo input file containing frequency responses of a 4-port network working in THz frequency range, which may be replaced by another touchstone file of interest.
- **Input_file_name:** The full name of input file including the file extention.
- **Num_poles:** The desired number of poles for vector fitting algorithm (i.e., vector fitting approximation order) in Full synthesis part of SROPEE.
- **VF_iter1:** Iteration number for calculating improved initial poles by fitting column sum through vector fitting algorithm in Full synthesis part of SROPEE.
- **VF_iter2:** Iteration number for calculating poles and residues through vector fitting algorithm in Full synthesis part of SROPEE.
- **weight_f:** Weighting type for vector fitting of function (f) in Full synthesis part of SROPEE.
- **weight_column_sum:** Weighting type for fitting of column sum in Full synthesis part of SROPEE.
- **VF_relax:** Enables/disables relaxed vector fitting in Full synthesis part of SROPEE.
- **VF_stable:** Enables/disables enforcing stability of fitted poles in Full synthesis part of SROPEE.
- **VF_asymp:** Fitting model options (fit with None, fit with D, fit with D and E) in Full synthesis part of SROPEE.
- **Passivity_Enforcement_Enable:** Enables/disables passivity assessment/enforcement in Full synthesis part of SROPEE.
- **SMP_iter_upper_limit:** Passivity enforcement iteration upper limit in Full synthesis part of SROPEE.
- **n:** The order of reducion for Bloack SAPOR algorithm in Passive MOR part of SROPEE.
- **in_port:** Port of input source for reverse MNA algorithm in Reduced synthesis part of SROPEE.
**Note:** The first input is located in folder: .\SROPEE\Input, and the remaining inputs are exported as input_variables.npy file in folder:. \SROPEE\Output.  



## Outputs:
- **input_variables.npy:** Input variables exported from SROPEE main program file (run_main.py) and defined by user. This file is used for transferring the user's input variables to different parts of SROPEE algorithm.
- **full_netlist.sp:** Equivalent full-order netlist file exported from Full synthesis part of SROPEE.
- **reduced_MNA_data.npz:** Reduced-order MNA matricres and data file exported from Passive MOR part of SROPEE.
- **reduced_netlist.sp:** Equivalent reduced-order netlist file exported from Reduced synthesis part of SROPEE.
- **reduced_voltages_ADS.cti:** Reduced-order nodal voltages obtained from simulation of equivalent reduced-order netlist in Keysight ADS software. This file is used for verification of Z-parameter (note: this file is generated with assumption of current source only in port 1, therefore, it will be used for plotting and comparing Z-parameters with inputs in port 1, (i.e, Z11, Z21, Z31, Z41). It may be replaced by another *.cti files generated with different assumptions).
- **figures:** Generated graphs for comparison and verification of obtained results in each part of SROPEE.


# Usage:
This program was designed to provide a reduced-order passivity-enforced equivalent circuit for the given S-parameter results of an optical interconnect.

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
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects," National Science Foundation, Award #1816542. Jun. 22, 2018.

[2] Guiana, Brian, GitHub. May 14, 2021. Available:https://github.com/bmguiana/OIDT

[3] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

[4] B. Gustavsen, "Improving the pole relocating properties of vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

[5] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method", IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.
     
[6] Simon De Ridder, GitHub. Feb 08, 2019. Accessed on: August 10, 2022, [Online].https://github.com/SimonDeRidder/pyVF

[7] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery, vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[8] Totorica, Nathan, GitHub, May 5, 2021. Accessed on: August 10, 2022, [Online]. Available:https://github.com/ntotorica/SMP_Passivity_Enforcement

[9] A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization through equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119

[10] V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from Frequency Response of Multi-port Networks", EMC EUR, pp. 79-84, 2016. [Online]. Available: ://WOS:000392194100012.

[11] R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent circuits of improper rational transfer functions with real poles and residues. IEEE Access, vol.8, page.125149:125162, 2020.

[12] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021, [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis

[13] http://scikit-rf.org

[14] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[15] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley and Sons, 2015.

[16] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou. SAPOR: second-order arnoldi method for passive order reduction of RCS circuits. In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.

[17] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou. Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

```
