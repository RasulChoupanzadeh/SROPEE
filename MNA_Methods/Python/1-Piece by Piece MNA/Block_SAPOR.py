
''' Block_SAPOR.py

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
from SOrth import SOrth

# Prepating input matrices
start = time.time()
D = 2*s0*Cm + Gm;
K = (s0**2)*Cm + s0*Gm + Gamma                          
B0 = s0*Bm
B1 = Bm
K_inv = np.linalg.inv(K)
x = np.hstack((-K_inv@D, K_inv))
y = np.hstack((-Cm, np.zeros(shape=(N,N))))
A = np.vstack((x, y))                          
Q0 = K_inv@B0;
P0 = B1

# Block SAPOR main algorithm (inputs: [A],[Q0],[P0],n,p     output: projection matrix [Q])
QP0 = np.vstack((Q0, P0))
k=int(n/p)    
QP = np.zeros(shape=(2*N,n), dtype='complex')
QP[:,0:p] = SOrth(QP0,p,N)
QP_hat = np.zeros(shape=(2*N,n), dtype='complex') 
for i in range(0,k-1):
    QP_hat[:, (i+1)*p:(i+2)*p] = A @ QP[:,i*p:(i+1)*p]
    for j in range(0,i+1):
        H = np.transpose(QP[0:N,j*p:(j+1)*p]) @ QP_hat[0:N,(i+1)*p:(i+2)*p]
        QP_hat[:,(i+1)*p:(i+2)*p] = QP_hat[:,(i+1)*p:(i+2)*p] - QP[:,j*p:(j+1)*p] @ H
    QP[:, (i+1)*p:(i+2)*p] = SOrth(QP_hat[:,(i+1)*p:(i+2)*p],p,N)
Q = QP[0:N,:]

# Calcualting reduced MNA matrices projected by [Q]
Cmr = np.transpose(Q) @ Cm @ Q
Gmr = np.transpose(Q) @ Gm @ Q
Gammar = np.transpose(Q) @ Gamma @ Q     
Bmr = np.transpose(Q) @ Bm
Lvr = np.transpose(Q) @ Lv 

# Calculating frequency response
V_reduced = np.zeros(shape=(q,len(freq)), dtype = 'complex')   
nt=0
for f in freq:
    s = 1j*2*np.pi*f*Coef
    V = np.transpose(Lvr) @ np.linalg.inv(s*Cmr + Gmr + Gammar/s) @ Bmr @ Jm
    V_reduced[:,nt]=V
    nt = nt+1

end = time.time()
