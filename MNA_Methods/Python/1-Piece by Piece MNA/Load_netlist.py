
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


# define function to call sub-circuits
def Load_sbck(sbck):
    Rabr = []; Lbr = []; Rbbr = []; Cbr = [];
    Num_branch = 0;
    j = sbck-1;
    for i in range(sbck_ident[j], sbck_ident[j+1]):
        
        if netlist[i][0:8]=='* Branch':
            Num_branch = Num_branch+1;
            
        if netlist[i][0:2]=='Ra':
            sentence = netlist[i];
            for t in sentence.split():
                try:
                    R = float(t)
                except ValueError:
                    pass
            Rabr.append(float(R))
            
        if netlist[i][0:2]=='Rb':
            sentence = netlist[i];
            for t in sentence.split():
                try:
                    R = float(t)
                except ValueError:
                    pass
            Rbbr.append(float(R))  
            
        if netlist[i][0:2]=='Lb':
            sentence = netlist[i];
            for t in sentence.split():
                try:
                    l = float(t)
                except ValueError:
                    pass
            Lbr.append(float(l))   
        if netlist[i][0:2]=='Cb':
            sentence = netlist[i];
            for t in sentence.split():
                try:
                    c = float(t)
                except ValueError:
                    pass
            Cbr.append(float(c))              
    Ra = np.array(Rabr);
    Rb = np.array(Rbbr);
    L = np.array(Lbr);
    C = np.array(Cbr);
    return Num_branch, Ra, L, Rb, C



