
''' general_MNA_builder.py

Author: Rasul Choupanzadeh
Date: 05/12/2022


This code is based on the concepts from [1-2].

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975

[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

'''



import numpy as np

def MNA_calculator(Ra, Rb, L, C, Num_branch, Num_port):
    N = 2*np.sum(Num_branch) + Num_port;
    M = np.sum(Num_branch);
    Ga = 1/Ra[:];
    Gb = 1/Rb[:];
    Gm = np.zeros(shape=(N,N));
    for i in range(0,Num_port):
        Gm[i,i] = 1;
        
    p = Num_port;       
    nt=0;
    for i in range (p,N,2):
        Gm[i,i] = Ga[nt];
        Gm[i+1,i+1] = Gb[nt];
        nt = nt+1;
    
    #Identify the number of branches connected to the ith port
    Num_br = Num_branch.copy();
    br_sbckix = np.zeros(shape=(p,1))
    t = p
    nt = 0
    while t != 0:
        br_sbckix[nt] = np.sum(Num_br[0:t])
        Num_br = Num_br[t:]
        t = t-1
        nt = nt+1
    
    # Correct the ith column of Gm (connected Ra (or Ga) to ith Node)    
    nt = p
    for j in range(0,p):
        for i in range(nt,nt+2*int(br_sbckix[j,0]),2):
            Gm[i,j] = -Gm[i,i]
        nt = nt + 2*int(br_sbckix[j,0])
        
        
    ## Correct the ith column of Gm (connected Rb (or Gb) to the ith Node)
    nt = p + 2*Num_branch[1]
    nk = 0
    t = p-1
    for k in range(1,p):
        for j in range(k,p):
            for i in range(nt,nt+2*Num_branch[j+nk],2):
                Gm[i+1,j] = -Gm[i+1,i+1]
            nt = nt+2*Num_branch[j+nk]
        nt = nt+2*Num_branch[j+nk+1]
        nk = nk+1
        t = t-1
    
    
    # Creating Cm matrix. Note: we treat Cm as same as Gb (or Rb) part of Gm.
    Cm = np.zeros(shape=(N,N))
    nt = 0
    for i in range(p,N,2):
        Cm[i+1,i+1] = C[nt]
        nt = nt+1
        
    ## Correct the ith column of Gm (connected Rb (or Gb) to the ith Node)
    nt = p + 2*Num_branch[1]
    nk = 0
    t = p-1
    for k in range(1,p):
        for j in range(k,p):
            for i in range(nt,nt+2*Num_branch[j+nk],2):
                Cm[i+1,j] = -Cm[i+1,i+1]
            nt = nt+2*Num_branch[j+nk]
        nt = nt+2*Num_branch[j+nk+1]
        nk = nk+1
        t = t-1    
        
        
    Lm = np.diag(L);
    
    Es = np.zeros(shape=(N,M));
    nt = 0;
    for i in range(p,N-1,2):
        Es[i,nt] = 1;
        Es[i+1,nt] = -1;
        nt = nt+1;
        
    SEt = np.linalg.inv(Lm) @ np.transpose(Es);
    Gamma = Es @ SEt;    
    
    Bm = np.zeros(shape=(N,p));
    nt = 0;
    for i in range (0,p):
        Bm[i,i] = 1
        
    Lv = np.identity(N)       
    return Gm, Cm, Gamma, Lv, Bm, M, N , SEt


