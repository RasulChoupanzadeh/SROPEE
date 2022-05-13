# Piece-by-Piece and General MNA Approaches for Model Order Reduction of Optical multiport Networks via SAPOR/Block SAPOR Algorithms
### Author: Rasul Choupanzadeh
### Date: 5/12/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, as part of the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1]. For the convenience of readers, the codes are written in both Python and MATLAB software.

# Overview
- There is a code that implements the MNA (Modified Nodal Analysis) approach descibed in [2-3], and constructs the MNA matrices due to the given SPICE netlist of multiport networks including canoniacal Ra-L-Rb-C and/or R-L circuit branches [4-5].

- It performs the MNA algorithm based on two proposed approaches:
    * **Piece-by-Piece MNA approach:** Divides the multiport network into smaller pieces (sub-circuits), then, performs MNA algorithm on each piece to obtain the Piece-by-Piece MNA system.
    * **General MNA approach:** Considers the entire network as a cohesive unit, and performs MNA algorithm on general network to obtain the General MNA system.

- It performs analysis on the full-order systems (system before order reduction) obtained by Piece-by-Piece and General MNA approaches, and calculates the Y-parameters (admittance parameters) of full-order netwrok in each approach.

- There are also codes that implement the SAPOR (Second-Order Arnoldi Method for Passive Order Reduction) and/or Block SAPOR algorithms [6-7], for Model Order Reduction (MOR) of multiport network, based directly on the presented algorithms in [6-7]. 

- They perform SAPOR and/or Block SAPOR techniques on the full-order systems obtained by Piece-by-Piece and General MNA approaches to reduce the orders of systems.

- It performs analysis on the reduced-order systems, which were obtained by applying SAPOR/Block SAPOR algorithms on the full-order Piece-by-Piece and General MNA matrices, and calculates the Y-parameters (admittance parameters) of the reduced-order netwrok in each approach.

- Finally, it provides the comparative results by calculating:
    * Full-order and reduced-order Y-parameters through Piece-by-Piece MNA approach
    * Full-order and reduced-order Y-parameters through General MNA approach
    * Execution time of solving full-order vs. reduced-order systems in Piece-by-Piece MNA approach
    * Execution time of solving full-order vs. reduced-order systems in General MNA approach
    * Root Mean Square Error (RMSE) of reduced-order and full-order networks in Piece-by-Piece MNA approach
    * Root Mean Square Error (RMSE) of reduced-order and full-order networks in General MNA approach
    * Oders of reduction in Piece-by-Piece and General MNA approaches


# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files:
**Note:** As with Python, MATLAB program files (.m) are named similarly to Python files (.py).
- **netlist.sp:** Input SPICE netlist file of a 4-port network operating in THz frequencies, could be replaced by other multiport network's netlist extracted from circuit simulators such as SPICE, Keysight ADS (Advanced Design System) or any other code/software.
- **Load_netlist.py:** Loads the SPICE netlist file, and then: 
    * Calculates the number of ports, sub-circuits, and branches
    * Stores the values of circuit elements (Ra, L, Rb, and C) in the corresponding arrays
- **MNA_calculator.py** and **general_MNA_builder.py:** Performs the MNA algorithm through Piece-by-Piece MNA and General MNA approaches, respectively.
- **original_system.py:** Performs analysis on the full-order system (original system) to calculate the nodal voltages. 
- **Block_SAPOR.py:** Calculates projection matrix using SAPOR/Block SAPOR algorithm, implements model order reduction, performs analysis of reduced-order system, and calculates the nodal voltages.    
- **SOrth.py:** A subfunction used for calculation of projection matrix.  
- **Run_Admittance_calculator.py:** It is the main program file, calculates the Y-parameters of full-order and reduced-order systems, RMSE, reduced system matrices, and execution times by calling and/or running all the previously mentioned program files.



# Full program
Run_Admittance_calculator.py in your desired approach (Piece-by-Piece MNA or General MNA approach). This code:
1. Runs Load_netlist. This is the script for loading netlist.sp file and calculating input parameters for MNA_calculator.py.
2. Runs MNA_calculator. This is the algorithm for calculating MNA matrices.
3. Runs original_system. This is the script for solving full-order system to find the nodal voltages.
4. Runs Block_SAPOR. This script performs SAPOR/Block SAPOR algorithm to reduce the order of system. 
5. Provides the reduced-order system matrices.
6. Calculates the RMSE between reduced-order and full-order systems.
7. Calculates the required time for solving full-order and reduced-order systems. 



# Run instructions
To run the algorithm for example of 4-port network operating in THz frequencies execute the Run_Admittance_calculator.py. To run the algorithm for your own network:
1. Replace your SPICE netlist file with existing netlist.sp.
2. Change the input parameters in Run_Admittance_calculator.py of your desired approach (Piece-by-Piece MNA or General MNA approach). The changeable input parameters are:
    * Coef: Normalization coefficient. It normalizes the order of operating frequency in KHz, MHz, GHz, THz, etc.
    * fmin: Minimum frequency (normalized by coefficient)
    * fmax: Maximum frequency (normalized by coefficient)
    * fstep: Frequency step (normalized by coefficient)
    * n0 and/or n_arr: This parameters exist only in Piece-by-Piece MNA approach. You can determine the reduction order for each sub-circuit either by using a fixed reduction order (n0) or by definig an array of reduction orders for sub-circuits (n_arr).
    * n: This parameter exists only in General MNA approach. It specifies the general order of reduction.
    * s0: Selected frequency expansion point (matching point). The projection matrix is created based on this maching point.
3. Run the program.


## Inputs:
- **Coef:** A scalar value for normalization coefficient
- **fmin:** A scalar value for notmalized minimum frequency point
- **fmax:** A scalar value for notmalized maximum frequency point
- **fstep:** A scalar value for notmalized frequency step
- **n0:** A scalar value for creating a fixed (identical) vector of reduction orders for sub-circuits in Piece-by-Piece MNA approach
- **n_arr:** Vector of reduction orders for sub-circuits in Piece-by-Piece MNA approach created either by defining n0 (fixed redcution orders), or by definig an array of reduction orders (non-identical order of reductions) 
    * Dimensions: number of sub-circuits (Num_sbck)
- **s0:** Selected frequency expansion point in operating frequency range
- **str_in:** An input string in the form of "**Ymn**", required to be entered by user, for comparing the desired Y-parameter. 
    * **m** & **n:** Integer values from 1 to number of ports (Num_port)

    
## Outputs:
- **Figure:** Plots the desired Y-parameter over frequency range:
    * Solid line: Y-parameter of full-order (original) system
    * Dotted line: Y-parameter of reduced-order system
- **RMSE:** Shows the Root Mean Square Error between the plotted Y-parameter of reduced-order and full-order systems.
- **t_full:** Execution time for solving full-order system.
- **t_reduced:** Execution time for solving reduced-order system.


# Usage:
This program was designed to provide a comparison of two proposed approaches for MNA analysis, Piece-by-Piece MNA and General MNA, for full-order (original) system and their efficiency to obtain the reduced-order system using SAPOR/Block SAPOR technique. The detailed comparison of results, advantages and disadvantages, and limitations and drawbacks are provided in a research article written by Rasul Choupanzadeh and Prof. Ata Zadehgol, which will be cited here once it has been published.


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
