# Passive Model Order Reduction (MOR) part of SROPEE algorithm
### Author: Rasul Choupanzadeh
### Date: 08/15/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This program uses General MNA approach for Model Order Reduction of Optical multiport Networks via Block SAPOR Algorithms. There is a code that implements the MNA (Modified Nodal Analysis) approach descibed in [2-3], and constructs the MNA matrices for a given SPICE netlist of multiport networks including canoniacal Ra-L-Rb-C and/or R-L circuit branches [4-5].

- It performs the MNA algorithm on the networks with the assumption of current input sources, which is a prerequisite condition for the reverse MNA part of SROPEE, based on the General MNA approach.
    * General MNA approach: Considers the entire network as a cohesive unit, and performs MNA algorithm on general network to obtain the General MNA system.
    * Note: A comparison of two methods (Piece-by-Piece and General MNA Approaches) is provided in ./SROPEE/MNA_Methods.

- It performs analysis on the full-order systems (system before order reduction) obtained by General MNA approaches, and calculates the Z-parameters of full-order netwrok.

- There is a code that implement the Block SAPOR (Block Second-Order Arnoldi Method for Passive Order Reduction) [6-7], for Model Order Reduction (MOR) of multiport network, based directly on the presented algorithms in [6-7]. 

- It performs Block SAPOR algorithm on the full-order system obtained by General MNA approaches to reduce the orders of systems.

- It performs analysis on the reduced-order systems, which were obtained by applying Block SAPOR algorithm on the full-order General MNA matrices, and calculates the Z-parameters of reduced-order netwrok.

- Finally, it provides the comparative results by calculating:
    * Full-order and reduced-order Z-parameters through General MNA approach
    * Execution time of solving full-order vs. reduced-order systems
    * Root Mean Square Error (RMSE) of reduced-order and full-order networks



# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files:
- **full_netlist.sp:** Input SPICE netlist file of optical interconnect extracted from Full synthesis part of SROPEE. (note: this file in located in the folder: .\SROPEE\Output, and may be replaced by other multiport network's netlist generated from Full synthesis part of SROPEE or circuit simulators such as SPICE and Keysight ADS).
- **Load_netlist.py:** Loads the SPICE netlist file, and then: 
    * Calculates the number of ports, sub-circuits, and branches.
    * Stores the values of circuit elements (Ra, L, Rb, and C, Rd) in corresponding arrays.
- **general_MNA_builder.py:** Performs the MNA algorithm through General MNA approach.
- **Block_SAPOR.py:** Implements model order reduction, calculates projection matrix using Block SAPOR algorithm, performs analysis of reduced-order system, and calculates reduced-order MNA matrices.    
- **SOrth.py:** A subfunction used for calculation of projection matrix.  
- **Run_passive_MOR.py:** This is the main program file for Passive MOR part of SROPEE. It calculates the Z-parameters of full-order and reduced-order systems, RMSE, full-order MNA matrices, reduced-order MNA matrices, and execution times, by executing the previously mentioned program files.
 


# Run instructions
This program is a part of SROPEE algorithm, and will be executed auromatically by executing the main program file of SROPEE (run_main.py)


## Inputs:
- **freq:** A vector of Frequency points (note: this vector is loaded from input file of SROPEE). The following variables will be obtained from this vector
    * fmin: Scalar value of minimum frequency point
    * fmax: Scalar value of maximum frequency point
    * fstep: Scalar value of frequency step
- **n:** Scalar value for order of redution (note: this input is loaded from input_variables.npy, which is located in the folder: .\SROPEE\Output)
- **s0:** Selected frequency expansion point in operating frequency range (note: this input variable be changed depending on user preference and frequency response of network)
- **Coef:** Scalar value for normalization of frequencies (note: this input variable is fixed in THz for optical interconnects of SROPEE algorithm, and may be changed depending on the operating frequency range)

    
## Outputs:
- **reduced_MNA_data.npz:** Reduced-order MNA matricres and data (note: this file is saved in the folder: .\SROPEE\Output)
- **RMSE:** The Root Mean Square Error between the reduced-order and full-order systems
- **t_full:** Execution time for solving full-order system.
- **t_reduced:** Execution time for solving reduced-order system.
- **Figure:** Z-parameter plots of original full-order vs. reduced-order networks shown in magnitude and phase (note: The graphs are saved as .pdf files in the folder: .\SROPEE\Output\figures\Passive_MOR).




# Usage:
This program was designed to implement the passive MOR on a given equivalent netlist of optical interconnects obtained from the first step of SROPEE. This step is considered as the second part of SROPEE algorithm.


# Software Version Information:
**Python 3.8.13**


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

[6] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou. SAPOR: second-order arnoldi method for passive order reduction of RCS circuits. In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pp. 74–79, 2004.

[7] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou. Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pp. 244–249, 2005.

```
