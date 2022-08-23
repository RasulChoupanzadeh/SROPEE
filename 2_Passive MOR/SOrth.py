
""" SOrth.py      => This script is a sub-function used in Block SAPOR algorithm for orthonormalization process.

Author: Rasul Choupanzadeh
Date: 07/03/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [3].

This code is based on the concepts from [1-2].

[1] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74_79, 2004.
    
[2] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244_249, 2005.
    
[3] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.
    
"""


 ## Input: QPm_hat,p,N           Output: qp


import numpy as np

def SOrth(QPm_hat,p,N):
    qp = np.zeros(shape=(2*N,p))
    R = np.zeros(shape=(p,p))
    for i in range(0,p):
        qp[:,i] = QPm_hat[:,i]
        for j in range(0,i):
            R[j,i] = np.transpose(qp[0:N,j]) @ qp[0:N,i]
            qp[:,i] = qp[:,i]-R[j,i]*qp[:,j]
        R[i,i] = np.linalg.norm(qp[0:N,i])
        if R[i,i] == 0:
            break
        qp[:,i] = qp[:,i]/R[i,i]
    return qp
    