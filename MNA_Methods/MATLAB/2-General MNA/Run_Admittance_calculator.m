%---------------------------------------------------------------------------------------------------------------------------------------%

% Run_Admittance_calculator.m
% 
% Author: Rasul Choupanzadeh
% Date: 05/12/2022
% 
% 
% This code uses MNA_calculator.m and Block_SAPOR.m, which are based on the concepts from [1-2] and [3-4], respectively.
% 
% [1] Chung-Wen Ho, A. Ruehli, and P. Brennan. The modified nodal approach to
%     network analysis. IEEE Transactions on Circuits and Systems, 22(6):504–509, 1975
% 
% [2] S. Grivet-Talocia and B. Gustavsen. Passive Macromodeling: Theory and Applications, 1st edition. John Wiley & Sons, 2015.
%
% [3] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
%     SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
%     In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.
% 
% [4] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
%     Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
%     RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

%---------------------------------------------------------------------------------------------------------------------------------------%

clear all
clc
% input parameters  
Coef=1e12;                       % Normalized coefficient for KHz               
fmin=175;                      
fmax=215;                      
fstep=0.1;                    
f=fmin:fstep:fmax;

n=200;                           % Order of reduction
s0=2*pi*194*Coef; 

% Build general MNA
run('general_MNA_builder.m')

prompt=sprintf('Please insert the desired atmittance as Ymn for m,n = 1,2,3 (e.g., Y11, Y12,...): ');
str1=input(prompt,'s');
in_port=str2double(str1(2));
out_port=str2double(str1(3));
Jm=zeros(p,1);
Jm(in_port)=1;

str2={'Y11','Y12','Y13','Y14';'Y21','Y22','Y23','Y24';'Y31','Y32','Y33','Y34';'Y41','Y42','Y43','Y44'};
Ymn = strcmp(str1,str2);

% Run original full order network from MNA
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
% Calculating current of inductances
V_node=V_orig(p+1:end,:);
I_ind=zeros(M,length(f));
nt=1;
for freq=fmin:fstep:fmax
    s=1i*2*pi*freq*Coef;
    I_ind(:,nt)=(SEt(:,p+1:N)*V_node(:,nt))./s;                         % Calculating currents through formula 3 of SAPOR paper
    nt=nt+1;
end 

I_sbck=zeros(Num_sbck,length(f));
nt=0;
for i=1:Num_sbck
    I_sbck(i,:)=sum(I_ind(nt+1:nt+Num_branch(i),:));
    nt=nt+Num_branch(i);
end

Y12=I_sbck(2,:);
Y13=I_sbck(3,:);
Y14=I_sbck(4,:);
Y11=I_sbck(1,:)+Y12+Y13+Y14;

Y23=I_sbck(6,:);
Y24=I_sbck(7,:);
Y22=I_sbck(5,:)-Y12+Y23+Y24;

Y34=I_sbck(9,:);
Y33=I_sbck(8,:)-Y13-Y23+Y34;

Y44=I_sbck(10,:)-Y14-Y24-Y34;

% Run Block SAPOR Algorithm for reduction
Num_port=p;
run('Block_SAPOR.m')
V_node_reduced=V_reduced(Num_port+1:end,:);
I_ind_reduced=zeros(M,length(f));
nt=1;
for freq=fmin:fstep:fmax
    s=1i*2*pi*freq*Coef;
    I_ind_reduced(:,nt)=(SEt(:,Num_port+1:N)*V_node_reduced(:,nt))./s;                         % Calculating currents through formula 3 of SAPOR paper
    nt=nt+1;
end 

I_sbck_reduced=zeros(Num_sbck,length(f));
nt=0;
for i=1:Num_sbck
    I_sbck_reduced(i,:)=sum(I_ind_reduced(nt+1:nt+Num_branch(i),:));
    nt=nt+Num_branch(i);
end

Y12_reduced=I_sbck_reduced(2,:);
Y13_reduced=I_sbck_reduced(3,:);
Y14_reduced=I_sbck_reduced(4,:);
Y11_reduced=I_sbck_reduced(1,:)+Y12_reduced+Y13_reduced+Y14_reduced;

Y23_reduced=I_sbck_reduced(6,:);
Y24_reduced=I_sbck_reduced(7,:);
Y22_reduced=I_sbck_reduced(5,:)-Y12_reduced+Y23_reduced+Y24_reduced;

Y34_reduced=I_sbck_reduced(9,:);
Y33_reduced=I_sbck_reduced(8,:)-Y13_reduced-Y23_reduced+Y34_reduced;

Y44_reduced=I_sbck_reduced(10,:)-Y14_reduced-Y24_reduced-Y34_reduced;

% Comparison of Original system (Full order) & Reduced system
if Ymn(1,1) == 1
    Y_full=Y11;
    Y_reduced=Y11_reduced;
elseif Ymn(1,2) == 1 || Ymn(2,1) == 1
    Y_full=Y12;
    Y_reduced=Y12_reduced;
elseif Ymn(1,3) == 1 || Ymn(3,1) == 1
    Y_full=Y13;
    Y_reduced=Y13_reduced;
elseif Ymn(1,4) == 1 || Ymn(4,1) == 1
    Y_full=Y14;
    Y_reduced=Y14_reduced;
elseif Ymn(2,2) == 1
    Y_full=Y22;
    Y_reduced=Y22_reduced;
elseif Ymn(2,3) == 1 || Ymn(3,2) == 1
    Y_full=Y23;
    Y_reduced=Y23_reduced;
elseif Ymn(2,4) == 1 || Ymn(4,2) == 1
    Y_full=Y24;
    Y_reduced=Y24_reduced;
elseif Ymn(3,3) == 1
    Y_full=Y33;
    Y_reduced=Y33_reduced;
elseif Ymn(3,4) == 1 || Ymn(4,3) == 1
    Y_full=Y34;
    Y_reduced=Y34_reduced;
elseif Ymn(4,4) == 1
    Y_full=Y44;
    Y_reduced=Y44_reduced;
end


semilogy(f,abs(Y_full),'b','LineWidth',2);
hold on
semilogy(f,abs(Y_reduced),'--r','LineWidth',2)

legend('Full order','Reduced order')
xlabel('f (THz)') 
ylabel(sprintf('| Y_{%s} |',str1(2:3))) 

RMSE = sqrt(mean((abs(Y_full)-abs(Y_reduced)).^2))
title(sprintf('Full order (N=504) vs Reduced order (n = %d) with RMSE = %d',n,RMSE))
grid on
