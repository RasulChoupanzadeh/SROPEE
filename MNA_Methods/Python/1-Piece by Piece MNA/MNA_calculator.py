
''' MNA_calculator.py

Author: Rasul Choupanzadeh
Date: 05/12/2022


This code is based on the concepts from [1-2].

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

'''



import numpy as np

def MNA_calculator(Ra, Rb, L, C, Num_branch):
    N = 2*Num_branch+1;
    M = Num_branch;
    Ga = 1/Ra[:];
    Gb = 1/Rb[:];
    Gm = np.zeros(shape=(N,N));
    Gm[0,0] = Ga.sum();
    nt=0;
    for i in range (1,N,2):
        Gm[0,i] = -Ga[nt];
        Gm[i,0] = -Ga[nt];
        Gm[i,i] = Ga[nt];
        Gm[i+1,i+1] = Gb[nt];
        nt = nt+1;
        
    Cm = np.zeros(shape=(N,N));
    nt = 0;
    for i in range (2,N,2):
        Cm[i,i] = C[nt];
        nt = nt+1;
        
    Lm = np.diag(L);
    
    Es = np.zeros(shape=(N,M));
    nt = 0;
    for i in range(1,N-1,2):
        Es[i,nt] = 1;
        Es[i+1,nt] = -1;
        nt = nt+1;
    
    # Voltage input corrections
    Es = np.delete(Es, 0, axis=0);      ## delete the first row of Es
    SEt = np.linalg.inv(Lm) @ np.transpose(Es);
    Gamma = Es @ SEt;
    
    Gm = np.delete(Gm, 0, axis=0);      ## delete the first row of Gm
    Cm = np.delete(Cm, 0, axis=0);      ## delete the first row of Cm
    Gm = np.delete(Gm, 0, axis=1);      ## delete the first column of Gm
    Cm = np.delete(Cm, 0, axis=1);      ## delete the first column of Cm    
    N = N-1;
    
    Bm = np.zeros(shape=(N,1));
    nt = 0;
    for i in range (0,N,2):
        Bm[i,0] = Ga[nt];
        nt = nt+1;
        
    Lv = np.zeros(shape=(N,M));
    nt = 0;
    for i in range(0,M):
        Lv[nt, i] = 1;
        nt = nt+2;
       
    return Gm, Cm, Gamma, Lv, Bm, M, N


