
""" original_system.py     => This script calculates frequency response of full-order (original) network using general MNA approach.

Author: Rasul Choupanzadeh
Date: 06/28/2022


This code is based on the concepts from [1-2].

[1] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pp. 74–79, 2004.
    
[2] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pp. 244–249, 2005.
    
"""


 ## Input: f, Cm, Gm, Gamma, Bm, Jm, Lv           Output:  V_orig
 

import numpy as np
import time

start = time.time()
V_orig = np.zeros(shape=(q,len(freq)), dtype = 'complex')
nt = 0
for f in freq:        
    s = 1j*2*np.pi*f*Coef
    V = np.transpose(Lv) @ (np.linalg.inv(s*Cm + Gm + Gamma/s) @ Bm) @ Jm
    V_orig[:,nt]=V[:,0]
    nt = nt+1

end = time.time()
