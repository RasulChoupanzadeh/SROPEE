
''' original_system.py

Author: Rasul Choupanzadeh
Date: 05/12/2022


This code is based on the concepts from [1-2].

[1] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
    SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
    In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.

[2] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
    Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
    RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

'''



# Calculate the original frequency response    
start = time.time()
V_orig = np.zeros(shape=(q,len(freq)), dtype = 'complex')
nt = 0
for f in freq:        
    s = 1j*2*np.pi*f*Coef
    V = np.transpose(Lv) @ (np.linalg.inv(s*Cm + Gm + Gamma/s) @ Bm) @ Jm
    V_orig[:,nt]=V[:,0]
    nt = nt+1

end = time.time()


