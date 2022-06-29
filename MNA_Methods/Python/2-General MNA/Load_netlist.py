
''' Load_netlist.py

Author: Rasul Choupanzadeh
Date: 05/12/2022

'''



import numpy as np

# Load and read the netlist line by line    
netlist = open('netlist.sp').readlines()                

# calculate the number of sub-circuits and ports
Num_sbck = 0;
sbck_ident = [];
for i in range(len(netlist)):
    if netlist[i][0]=='X':
        Num_sbck = Num_sbck+1
        Num_port = int(netlist[i][2])
    if netlist[i][0:4]=='.sub':
        sbck_ident.append(i)        
del sbck_ident[0]                       ## to delete a row from list==> del matrix[i,:]          to delete a row from a numpy array==> np.delete(array name, row/column number, axis=0/1 for row/column)
sbck_ident.append(len(netlist))


Num_branch =[];
for j in range(0,Num_sbck):
    Num_b = 0;
    for i in range(sbck_ident[j], sbck_ident[j+1]):
        if netlist[i][0:8]=='* Branch':
            Num_b = Num_b+1;
    Num_branch.append(Num_b)
    
Ra = []; 
L = []; 
Rb = []; 
C = []; 
for i in range(0,len(netlist)):
    if netlist[i][0:2]=='Ra':
        sentence = netlist[i];
        for t in sentence.split():
            try:
                R = float(t)
            except ValueError:
                pass
        Ra.append(float(R))
        
    if netlist[i][0:2]=='Rb':
        sentence = netlist[i];
        for t in sentence.split():
            try:
                R = float(t)
            except ValueError:
                pass
        Rb.append(float(R))  
        
    if netlist[i][0:2]=='Lb':
        sentence = netlist[i];
        for t in sentence.split():
            try:
                l = float(t)
            except ValueError:
                pass
        L.append(float(l)) 
        
    if netlist[i][0:2]=='Cb':
        sentence = netlist[i];
        for t in sentence.split():
            try:
                c = float(t)
            except ValueError:
                pass
        C.append(float(c))              


Num_branch = np.array(Num_branch);
Ra = np.array(Ra);
Rb = np.array(Rb);
L = np.array(L);
C = np.array(C);


