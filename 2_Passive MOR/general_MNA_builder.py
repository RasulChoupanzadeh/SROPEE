
""" general_MNA_builder.py

Author: Rasul Choupanzadeh
Date: 08/23/2022

# Acknowledgement
This project is completed as part of research conducted with my major professor and advisor, Prof. Ata Zadehgol, in the Applied Computational Electromagnetics and Signal/Power Integrity (ACEM-SPI) Lab while working toward the Ph.D. in Electrical Engineering at the University of Idaho, Moscow, Idaho, USA. 
This project was supported by a research grant from the National Science Foundation, under the NSF Award #1816542 [3].

This code is based on the concepts from [1-2].

[1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
    network analysis. IEEE Transactions on Circuits and Systems, 22(6):504_509, 1975.
    
[2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.

[3] A. Zadehgol, "SHF: SMALL: A Novel Algorithm for Automated Synthesis of Passive, Causal, and Stable Models for Optical Interconnects", National Science Foundation, Award #1816542. Jun. 22, 2018.

"""


## Input: Ra, Rb, L, C, Rd, Num_branch, Num_port           Output: Gm, Cm, Gamma, Lv, Bm, M, N , SEt


import numpy as np

def MNA_calculator(Ra, Rb, L, C, Rd, Num_branch, Num_port):
    if Rd.shape != 0:
        Num_branch = Num_branch - 1
        
    N = 2*np.sum(Num_branch) + Num_port;
    M = np.sum(Num_branch);       
    Ga = 1/Ra[:];
    Gb = 1/Rb[:];
    Gd = 1/Rd[:];
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
        
    # Correct the ith column of Gm (connected Rb (or Gb) to the ith Node)
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
        
    # Correct the ith column of Gm (connected Rb (or Gb) to the ith Node)
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
    
    
    
    # If the input sources are current source:
    #  - Es, Bm, Lm, Lv will remain unchanged. We just need to change the first p equations related to the input ports. Since the ports are connected only to R (Ra and Rb) and C elements, only the first p equations of Cm and Gm must be changed.
    #  - For this purpose, we must symmetrize Cm and Gm, and change the diagonal values of first p equations into the -sum of connected elements of that row (equation).
    #  - To symmetrize Cm,Gm ( which are lower triangular matrices), we can add the transpose of the lower triangular matrix with elements below the main diagonal to original matrix
    
    # Symmetrize Cm,Gm
    Gm2 = np.tril(Gm,-1);                   # Creates lower triangular of Gm with elements below the main diagonal
    Gm = Gm + np.transpose(Gm2);            # To calculate transpose of A, you can use both np.transpose(A) or A.T 
    
    Cm2 = np.tril(Cm,-1);                   
    Cm = Cm + np.transpose(Cm2);     
    
    # Change the diagonal values of first p equations
    for i in range (0,p):
        Gm[i,i] = 0
        Gm[i,i] = -np.sum(Gm[i,:])
        Cm[i,i] = 0
        Cm[i,i] = -np.sum(Cm[i,:])        
        
    
    ## Rd does not add new unknown vlotage, or auxiliary current==> N, and M will not be changed ==> only some values will be added to Gm[0:p,0:p].
    # Reshape Gd to p*p
    def create_upper_matrix(values, size):
        upper = np.zeros((size, size),dtype='complex')
        upper[np.triu_indices(size, 0)] = values
        return(upper)

    u = create_upper_matrix(Gd, Num_port).real
    l = np.triu(u,1).T
    Gd = u + l    
    
    # Calculated additional values Gmd, which shows the related equations between V1,...,Vp
    Gmd = np.zeros(shape=Gd.shape)
    for i in range (0,p):
        for j in range (0,p):
            if i==j:
                Gmd[i,j] = np.sum(Gd[i,:])  # or =np.sum(Gd,0) 
            else:
                Gmd[i,j] = -Gd[i,j]
    

    # Modify Gm by adding Gmd 
    for i in range (0,p):
        for j in range (0,p):
            Gm[i,j] = Gm[i,j] + Gmd[i,j]
             
        
    # To check the symmetry
    #if Gm.any() == Gm.T.any():
        #print('Gm is symmetric')
    #else:
        #print('Gm is not symmetric')
    
    return Gm, Cm, Gamma, Lv, Bm, M, N , SEt


