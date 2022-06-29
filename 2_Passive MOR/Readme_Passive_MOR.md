# Passive Model Order Reduction (MOR) of SROPEE
### Author: Rasul Choupanzadeh
### Date: 06/28/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program uses General MNA approach for Model Order Reduction of Optical multiport Networks via Block SAPOR Algorithms. There is a code that implements the MNA (Modified Nodal Analysis) approach descibed in [2-3], and constructs the MNA matrices due to the given SPICE netlist of multiport networks including canoniacal Ra-L-Rb-C and/or R-L circuit branches [4-5].

- It performs the MNA algorithm on the networks with current inputs, which is a prerequisite condition for reverse MNA part of SROPEE, based on the General MNA approach.
    * **General MNA approach:** Considers the entire network as a cohesive unit, and performs MNA algorithm on general network to obtain the General MNA system.

- It performs analysis on the full-order systems (system before order reduction) obtained by General MNA approaches, and calculates the Z-parameters of full-order netwrok.

- There are also codes that implement the Block SAPOR (Block Second-Order Arnoldi Method for Passive Order Reduction) [6-7], for Model Order Reduction (MOR) of multiport network, based directly on the presented algorithms in [6-7]. 

- It performs Block SAPOR algorithm on the full-order system obtained by General MNA approaches to reduce the orders of systems.

- It performs analysis on the reduced-order systems, which were obtained by applying Block SAPOR algorithm on the full-order General MNA matrices, and calculates the Z-parameters of reduced-order netwrok.

- Finally, it provides the comparative results by calculating:
    * Full-order and reduced-order Z-parameters through General MNA approach
    * Execution time of solving full-order vs. reduced-order systems in General MNA approach
    * Root Mean Square Error (RMSE) of reduced-order and full-order networks in General MNA approach
    * Oders of reduction in General MNA approach


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files:
- **full_netlist.sp:** Input SPICE netlist file of a 4-port network operating in THz frequencies, could be replaced by other multiport network's netlist extracted from Full synthesis part of SROPEE, or circuit simulators such as SPICE, Keysight ADS (Advanced Design System).
- **Load_netlist.py:** Loads the SPICE netlist file, and then: 
    * Calculates the number of ports, sub-circuits, and branches
    * Stores the values of circuit elements (Ra, L, Rb, and C) in the corresponding arrays
- **general_MNA_builder.py:** Performs the MNA algorithm through General MNA approach.
- **original_system.py:** Performs analysis on the full-order system (original system) to calculate the nodal voltages. 
- **Block_SAPOR.py:** Calculates projection matrix using Block SAPOR algorithm, implements model order reduction, performs analysis of reduced-order system, and calculates the nodal voltages.    
- **SOrth.py:** A subfunction used for calculation of projection matrix.  
- **Z_param_full.s4p:** Z-parameter results of full-order network simulated in ADS software.  
- **Load_Z_ADS_full.py:** Loads the results of full-order network simulated in ADS software into python for verification of the results of full and reduced order systems 
- **Run_main_passive_MOR.py:** It is the main program file, calculates the Z-parameters of full-order and reduced-order systems, RMSE, reduced system matrices, and execution times by calling and/or running all the previously mentioned program files.
 


# Full program
Execute the Run_main_passive_MOR.py script. This code:
1. Runs Load_netlist. This is the script for loading netlist.sp file and calculating input parameters for MNA_calculator.py.
2. Runs MNA_calculator. This is the algorithm for calculating MNA matrices.
3. Runs original_system. This is the script for solving full-order system to find the nodal voltages.
4. Runs Block_SAPOR. This script performs SAPOR/Block SAPOR algorithm to reduce the order of system. 
5. Provides the reduced-order system matrices.
6. Calculates the RMSE between reduced-order and full-order systems.
7. Calculates the required time for solving full-order and reduced-order systems. 



# Run instructions
To run the algorithm on the example of 4-port network operating in THz frequencies, execute the Run_Admittance_calculator.py. To run the algorithm for your own network:
1. Replace your SPICE netlist file with existing full_netlist.sp.
2. Change the input parameters in Run_main_passive_MOR.py. The changeable input parameters are:
    * Coef: Normalization coefficient. It normalizes the order of operating frequency in KHz, MHz, GHz, THz, etc.
    * fmin: Minimum frequency (normalized by coefficient)
    * fmax: Maximum frequency (normalized by coefficient)
    * fstep: Frequency step (normalized by coefficient)
    * n: This parameter specifies the general order of reduction.
    * s0: Selected frequency expansion point (matching point). The projection matrix is created based on this maching point.
3. Execute the Run_main_passive_MOR.py program.


## Inputs:
- **Coef:** A scalar value for normalization coefficient
- **fmin:** A scalar value for notmalized minimum frequency point
- **fmax:** A scalar value for notmalized maximum frequency point
- **fstep:** A scalar value for notmalized frequency step
- **n:** A scalar value for order of redution
- **s0:** Selected frequency expansion point in operating frequency range
- **str_in:** An input string in the form of "**Zmn**", required to be entered by user, for comparing the desired Z-parameter. 
    * **m** & **n:** Integer values from 1 to number of ports (Num_port)

    
## Outputs:
- **Figure:** Plots the desired Z-parameter over frequency range:
    * Solid line: Z-parameter of full-order (original) system
    * Dotted line: Z-parameter of reduced-order system
- **RMSE:** Shows the Root Mean Square Error between the plotted Z-parameter of reduced-order and full-order systems.
- **t_full:** Execution time for solving full-order system.
- **t_reduced:** Execution time for solving reduced-order system.
- **reduced_MNA_data.npz:** Reduced MNA matricres and data.


# Usage:
This program was designed to implement the passive MOR on a given equivalent netlist of optical interconnects obtained from the first step of SROPEE. This step is considered as the middle step (second step) of SROPEE.


# Software Version Information:
**Python 3.8.13**

**MATLAB R2014a**

Libraries used in Python:
   * pip		21.2.2
   * scikit-learn	1.0.2
   * numpy		1.21.5
   * matplotlib	        3.5.1




# References:
```
[1] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects," National Science Foundation, Award #1816542. Jun. 22, 2018.

[2] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[3] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley and Sons, 2015.

[4] A. Zadehgol. A semi-analytic and cellular approach to rational system characterization through equivalent circuits. International Journal of Numerical Modelling: Electronic Networks, Devices and Fields, 29(4):637–652, 2016.

[5] R. Choupanzadeh and A. Zadehgol. Stability, causality, and passivity analysis of canonical equivalent circuits of improper rational transfer functions with real poles and residues. IEEE Access, 8:125149–125162, 2020.

[6] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou. SAPOR: second-order arnoldi method for passive order reduction of RCS circuits. In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.

[7] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou. Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

```
