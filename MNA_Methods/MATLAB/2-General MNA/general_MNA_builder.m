%---------------------------------------------------------------------------------------------------------------------------------------%

% general_MNA_builder.m 
% 
% Author: Rasul Choupanzadeh
% Date: 05/12/2022
% 
% 
% This code is based on the concepts from [1-2].
% 
% [1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
%     network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975
% 
% [2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.
%
%
% Note: This code includes both the equivalent codes for Load_netlist.py and general_MNA_builder.py
%
%---------------------------------------------------------------------------------------------------------------------------------------%


% Load netlist
global netlist ;
file='netlist_passive.sp';
fileID = fopen(file);
netlist = textscan(fileID,'%s %s %s %s %s %s %s %s %s');
fclose(fileID);
clear fileID file ans;

%Split the sub-circuits
subckt_identifier=[];
nt = 1;
for i=1:length(netlist{1})
    s = netlist{1}{i};
    st= '.ends';
    tf = strcmp(s,st);                          % Compares the strings with '.ends'
    if tf==1
        subckt_identifier(nt)=i+1;              % Stores the row numbers for the begining of each sub_circuit; which is occures after each '.ends'
        nt=nt+1;
    end  
end
strg = netlist{1}{subckt_identifier(1)-2};      % Last sbck's name
p = str2double(strg(3));                        % Number of ports

Num_sbck=length(subckt_identifier)-1;           % Number of subcircuits 
Num_branch=zeros(Num_sbck,1);                   % Vector containing number of branches for each subcircuit                       
for i=1:Num_sbck
    Num_branch(i)=0;
    for j=subckt_identifier(i):subckt_identifier(i+1)-1
        s = netlist{2}{j};
        st= 'Branch';
        tf = strcmp(s,st);                      
        if tf==1
            Num_branch(i)=Num_branch(i)+1;              
        end  
    end
end

% Total number of voltage nodes (N)
N0=0;                      
for i=1:Num_sbck
    N0=N0+2*Num_branch(i);
end
N=N0+p;

% Total number of inductance currents (M)
M=sum(Num_branch); 

% Flattened vectors Ra,L,Rb,and C, containing pertinent elements of all branches of all subcircuits. 
Ra=zeros(1,M);
L=zeros(1,M);
Rb=zeros(1,M);
C=zeros(1,M);

nt1=1;
nt2=1;
nt3=1;
nt4=1;
for i=1:length(netlist{1})
    s = netlist{1}{i};
    switch(s(1))
        case{'R'}
            switch(s(2))
                case{'a'}
                    Ra(nt1)= str2double(netlist{4}{i});
                    nt1=nt1+1;
                case{'b'}
                    Rb(nt2)= str2double(netlist{4}{i});
                    nt2=nt2+1;
            end
        case{'L'}
            L(nt3)= str2double(netlist{4}{i});
            nt3=nt3+1;
        case{'C'}
            C(nt4)= str2double(netlist{4}{i});
            nt4=nt4+1;
    end
end
% Ra=abs(Ra);
% Rb=abs(Rb);
% C=abs(C);
% L=abs(L);

Ga=1./Ra;
Gb=1./Rb;

Gm=zeros(N,N);
for i=1:p
    Gm(i,i)=1;
end

nt=1;
for i=p+1:2:N
    Gm(i,i)=Ga(nt);
    Gm(i+1,i+1)=Gb(nt);
    nt=nt+1;
end

Num_br=Num_branch;
br_sbckix=zeros(1,p);
t=p;
nt=1;
while t~=0
    br_sbckix(nt)=sum(Num_br(1:t));               
    Num_br=Num_br(t+1:end);
    t=t-1;
    nt=nt+1;
end

nt=p;
for j=1:p
    for i=nt+1:2:nt+2*br_sbckix(j)
        Gm(i,j)=-Gm(i,i);                       % Corrects the ith column of Gm (connected Ra (or Ga) to ith Node)
    end
    nt=nt+2*br_sbckix(j);
end


nt=p+2*Num_branch(1);
nk=0;
t=p-1;
for k=2:p
    for j=k:p
        for i=nt+1:2:nt+2*Num_branch(j+nk)
            Gm(i+1,j)=-Gm(i+1,i+1);                 % Corrects the ith column of Gm (connected Rb (or Gb) to the ith Node)
        end
        nt=nt+2*Num_branch(j+nk);
    end
    nt=nt+2*Num_branch(j+nk+1);
    nk=nk+t;
    t=t-1;
end

% Creating Cm matrix. Note: we treat Cm as same as Gb (or Rb) part of Gm.
Cm=zeros(N,N);
nt=1;
for i=p+1:2:N
    Cm(i+1,i+1)=C(nt);
    nt=nt+1;
end

nt=p+2*Num_branch(1);
nk=0;
t=p-1;
for k=2:p
    for j=k:p
        for i=nt+1:2:nt+2*Num_branch(j+nk)
            Cm(i+1,j)=-Cm(i+1,i+1);                 % Corrects the ith column of Cm (connected C to the ith Node)
        end
        nt=nt+2*Num_branch(j+nk);
    end
    nt=nt+2*Num_branch(j+nk+1);
    nk=nk+t;
    t=t-1;
end

Lm=diag(L);
Es=zeros(N,M);
nt=1;
for i=p+1:2:N-1
    Es(i,nt)=1;
    Es(i+1,nt)=-1;
    nt=nt+1;
end
SEt=Lm\transpose(Es);
Gamma=Es*SEt;

Bm=zeros(N,p);
nt=1;
for i=1:p
    Bm(i,i)=1;
    nt=nt+1;
end

Lv=eye(N);

