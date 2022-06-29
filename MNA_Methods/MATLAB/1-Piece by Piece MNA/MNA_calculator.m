%---------------------------------------------------------------------------------------------------------------------------------------%

% MNA_calculator.m
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

%---------------------------------------------------------------------------------------------------------------------------------------%


function [Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch)
N=2*Num_branch+1;                      % Number of voltage nodes
M=Num_branch;                          % Number of inductance currents

Ga=1./Ra;
Gb=1./Rb;
Gm=zeros(N,N);
Gm(1,1)=sum(Ga);
nt=1;
for i=2:2:N
    Gm(1,i)=-Ga(nt);
    Gm(i,1)=-Ga(nt);
    Gm(i,i)=Ga(nt);
    Gm(i+1,i+1)=Gb(nt);
    nt=nt+1;
end

Cm=zeros(N,N);
nt=1;
for i=3:2:N
    Cm(i,i)=C(nt);
    nt=nt+1;
end

Lm=diag(L);

Es=zeros(N,M);
nt=1;
for i=2:2:N-1
    Es(i,nt)=1;
    Es(i+1,nt)=-1;
    nt=nt+1;
end
%-------------------------------------Voltage input-----------------------------
Es(1,:)=[];
SEt=Lm\transpose(Es);
Gamma=Es*SEt;

Gm(1,:)=[];
Gm(:,1)=[];
Cm(1,:)=[];
Cm(:,1)=[];
N=N-1;

Bm=zeros(N,1);
nt=1;
for i=1:2:N
    Bm(i,1)=1/Ra(nt);
    nt=nt+1;
end

Lv=zeros(N,M);
nt=1;
for j=1:M
    Lv(nt,j)=1;
    nt=nt+2;
end

end


