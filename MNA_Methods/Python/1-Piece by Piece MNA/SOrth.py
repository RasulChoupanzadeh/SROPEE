
''' SOrth.py

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



import numpy as np

def SOrth(QPm_hat,p,N):
    qp = np.zeros(shape=(2*N,p), dtype='complex')
    R = np.zeros(shape=(p,p), dtype='complex')
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
    

