# S-parameter to Reduced-Order Passivity-Enforced Equivalent Circuit (SROPEE)
### Author: Rasul Choupanzadeh
### Date: 06/28/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [1].

# Overview
- This algorithm has 3 main parts (3 main folder) in the following orders:
    1. Full synthesis
    2. Passive Model Order Reduction (Passive MOR)
    3. Reduced synthesis

- The Full synthesis part implements circuit synthesis on a given S-parameter file by performing:
    * Vector Fitting algorithm to obtain fitted poles/residues
    * Calculating circuit parameters from obtained poles/residues
    * Generating full SPICE netlist from obtained circuit parameters 

- The Passvie MOR part implements a model order reduction algorithm (Block SAPOR) on a given SPICE netlist file by performing:
    * General MNA method for calculating full MNA matrices
    * Block SAPOR algorithm for calculating reduced MNA matrices from full MNA matrices

- The Reduced synthesis part implements circuit synthesis on the given reduced MNA matrices by performing:
    * Reverse MNA algorithm to obtain reduced circuit parameters from reduced MNA matrices
    * Generating reduced SPICE netlist from obtained reduced circuit parameters 

- Finally, it provides an equivalent reduced netlist, which can be imported to Keysight ADS software for simulation and verification of obtained results.



# Licensing
Licensed under GNU GPL v.3.
 

# Files:

## Primary Program Files/Folders:
- **Full synthesis folder:** Implements full circuit synthesis.
    * Inputs: S_par.sxp
    * Outputs: full_netlist.sp 
- **Passive MOR:** Implements model order reduction algorithm.
    * Inputs: full_netlist.sp
    * Outputs: reduced_MNA_data.npz
- **Reduced synthesis folder:** Implements reduced circuit synthesis.
    * Inputs: reduced_MNA_data.npz
    * Outputs: reduced_netlist.sp 
- **Flow chart.pdf:** Provides the general flow chart of SROPEE algorithm.
- **Instruction.pdf:** Provides fully detailed steps for executing SROPEE algorithm. 


# Run instructions
To see the fully detailed instructions, please review the Instruction.pdf file.


# Usage:
This program was designed to provide an equivalent reduced netlist for a given S-parameter results of an optical interconnect.

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

[12] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[13] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley and Sons, 2015.

[14] A. Zadehgol. A semi-analytic and cellular approach to rational system characterization through equivalent circuits. International Journal of Numerical Modelling: Electronic Networks, Devices and Fields, 29(4):637–652, 2016.

[15] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou. SAPOR: second-order arnoldi method for passive order reduction of RCS circuits. In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.

[16] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou. Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

```

