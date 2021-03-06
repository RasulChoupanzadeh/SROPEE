# S-parameter to Reduced-Order Passivity-Enforced Equivalent Circuit (SROPEE)
### Author: Rasul Choupanzadeh
### Date: 07/03/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This algorithm has three main parts in three different folders in the following order:
    1. Full synthesis
    2. Passive Model Order Reduction (Passive MOR)
    3. Reduced synthesis

- The Full synthesis part implements a full-order circuit synthesis on a given S-parameter file (input file) by performing:
    * Vector Fitting algorithm to obtain fitted poles/residues
    * Calculating full-order circuit parameters from poles/residues
    * Generating full-order SPICE netlist from obtained circuit parameters 

- The Passvie MOR part implements Block SAPOR as the model order reduction technique on the full-order SPICE netlist file exported from Full synthesis part of SROPEE. It performs:
    * General MNA method for calculating full-order MNA matrices
    * Block SAPOR algorithm to calculate reduced-order MNA matrices

- The Reduced synthesis part implements circuit synthesis on the reduced-order MNA matrices exported from Passive MOR part. It performs:
    * Reverse MNA algorithm to obtain reduced-order circuit parameters from reduced-order MNA matrices
    * Generating reduced-order SPICE netlist from obtained reduced-order circuit parameters 

- Finally, it provides an equivalent reduced-order netlist, which may be imported to Keysight ADS software for further simulations and verification of obtained results.



# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Full synthesis folder:** Implements full-order circuit synthesis to generate equivalent full-order netlist.
    * Inputs: te_rough_example.s4p 
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
- **te_rough_example.s4p:** A demo input file of an optical interconnect frequency response, which may be replaced by another touchstone file of interest.
- **Input_file_name:** The full name of input file including the file extention.
- **Num_pole_pairs:** The number of pairs of poles for vector fitting algorithm in Full synthesis part of SROPEE.
- **n:** The order of reducion for Bloack SAPOR algorithm in Passive MOR part of SROPEE.
- **in_port:** Port of input source for reverse MNA algorithm in Reduced synthesis part of SROPEE.
**Note:** The first input is located in the folder: .\SROPEE\Input, and the remaining inputs are exported as input_variables.npy file in the folder:. \SROPEE\Output.  



## Outputs:
- **input_variables.npy:** Input variables exported from SROPEE main program file (run_main.py) and defined by user. This file is used for transferring the user's input variables to various parts of SROPEE algorithm.
- **full_netlist.sp:** Equivalent full-order netlist file exported from Full synthesis part of SROPEE.
- **reduced_MNA_data.npz:** Reduced-order MNA matricres and data file exported from Passive MOR part of SROPEE.
- **reduced_netlist.sp:** Equivalent reduced-order netlist file exported from Reduced synthesis part of SROPEE.
- **reduced_voltages_ADS.cti:** Reduced-order nodal voltages obtained from simulation of equivalent reduced-order netlist in Keysight ADS software. This file is used for verification of Z-parameter (note: this file is generated with the assumption of current source only in port 1, therefore, it will be used for plotting and comparing Z-parameters with inputs in port 1, (i.e, Z11, Z21, Z31, Z41). It may be replaced by another *.cti files generated with different assumptions).
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

[12] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to network analysis. IEEE Transactions on Circuits and Systems, 22(6):504???509, 1975

[13] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley and Sons, 2015.

[14] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou. SAPOR: second-order arnoldi method for passive order reduction of RCS circuits. In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74???79, 2004.

[15] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou. Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244???249, 2005.

```

