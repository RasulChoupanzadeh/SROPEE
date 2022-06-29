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
Coef=1e12;                       % Normalized coefficient for THz (changeable to MHz, GHz, etc)                  
fmin=175;                      
fmax=215;                      
fstep=0.1;                    
f=fmin:fstep:fmax;
Jm=[1];                         % Value of input source (in frequency domain)


% ---------------------------------------------Piece-wise calculation of output (currnet)-------------------------------------------
%sbck=1
n=35;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(1);                             % Loads sbck 1 (subcircuit 1) from netlist    
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 1
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig11=0;
for i=1:length(Ra)
    I_orig11=I_orig11+(V_orig(i,:)-Jm)/Ra(i);                        % Total current of voltage source for original sbck 1
end
I_red11=0;
for i=1:length(Ra)
    I_red11=I_red11+(V_reduced(i,:)-Jm)/Ra(i);                       % Total current of voltage source for reduced order sbck 1 
end 

%sbck=2
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(2);                             % Loads sbck 2 (subcircuit 2) from netlist                        
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 2 
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig12=0;
for i=1:length(Ra)
    I_orig12=I_orig12+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 2                
end
I_red12=0;
for i=1:length(Ra)
    I_red12=I_red12+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 2
end 

%sbck=3
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(3);                             % Loads sbck 3 (subcircuit 3) from netlist               
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 3 
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig13=0;
for i=1:length(Ra)
    I_orig13=I_orig13+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 3
end
I_red13=0;
for i=1:length(Ra)
    I_red13=I_red13+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 3
end

%sbck=4
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(4);                             % Loads sbck 4 (subcircuit 4) from netlist          
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 4
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig14=0;
for i=1:length(Ra)
    I_orig14=I_orig14+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 4
end
I_red14=0;
for i=1:length(Ra)
    I_red14=I_red14+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 4
end 


%sbck=5
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(5);                             % Loads sbck 5 (subcircuit 5) from netlist          
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 5
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig22=0;
for i=1:length(Ra)
    I_orig22=I_orig22+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 5
end
I_red22=0;
for i=1:length(Ra)
    I_red22=I_red22+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 5
end 

%sbck=6
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(6);                             % Loads sbck 6 (subcircuit 6) from netlist                      
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 6
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig23=0;
for i=1:length(Ra)
    I_orig23=I_orig23+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 6
end
I_red23=0;
for i=1:length(Ra)
    I_red23=I_red23+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 6
end 

%sbck=7
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(7);                             % Loads sbck 7 (subcircuit 7) from netlist                      
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 7
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig24=0;
for i=1:length(Ra)
    I_orig24=I_orig24+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 7
end
I_red24=0;
for i=1:length(Ra)
    I_red24=I_red24+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 7
end

%sbck=8
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(8);                             % Loads sbck 8 (subcircuit 8) from netlist            
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 8
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig33=0;
for i=1:length(Ra)
    I_orig33=I_orig33+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 8
end
I_red33=0;
for i=1:length(Ra)
    I_red33=I_red33+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 8
end

%sbck=9
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(9);                             % Loads sbck 9 (subcircuit 9) from netlist            
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 9
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig34=0;
for i=1:length(Ra)
    I_orig34=I_orig34+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 9
end
I_red34=0;
for i=1:length(Ra)
    I_red34=I_red34+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 9
end

%sbck=10
n=30;
s0=2*pi*194*Coef; 
[Ra,Rb,L,C,Num_branch]=Load_netlist(10);                             % Loads sbck 10 (subcircuit 10) from netlist            
[Gm,Cm,Gamma,Lv,Bm,M,N]=MNA_calculator(Ra,Rb,L,C,Num_branch);       % Calculates MNA matrices for sbck 10
p=size(Bm,2);
q=size(Lv,2);
tic
run('original_system.m')
toc
run('Block_SAPOR.m')
I_orig44=0;
for i=1:length(Ra)
    I_orig44=I_orig44+(V_orig(i,:)-1)/Ra(i);                        % Total current of voltage source for original sbck 10
end
I_red44=0;
for i=1:length(Ra)
    I_red44=I_red44+(V_reduced(i,:)-1)/Ra(i);                       % Total current of voltage source for reduced order sbck 10
end


% ---------------------------------------------Y-parameter calculation-------------------------------------------
Y11_original=(I_orig11)+(I_orig12)+(I_orig13)+(I_orig14);
Y11_reduced=(I_red11)+(I_red12)+(I_red13)+(I_red14);

Y21_original=(I_orig12);
Y21_reduced=(I_red12);

Y31_original=(I_orig13);
Y31_reduced=(I_red13);

Y41_original=(I_orig14);
Y41_reduced=(I_red14);

Y22_original=(I_orig12)+(I_orig22)+(I_orig23)+(I_orig24);
Y22_reduced=(I_red12)+(I_red22)+(I_red23)+(I_red24);

Y32_original=(I_orig23);
Y32_reduced=(I_red23);

Y42_original=(I_orig24);
Y42_reduced=(I_red24);

Y33_original=(I_orig13)+(I_orig23)+(I_orig33)+(I_orig34);
Y33_reduced=(I_red13)+(I_red23)+(I_red33)+(I_red34);

Y43_original=(I_orig34);
Y43_reduced=(I_red34);

Y44_original=(I_orig14)+(I_orig24)+(I_orig34)+(I_orig44);
Y44_reduced=(I_red14)+(I_red24)+(I_red34)+(I_red44);

% ---------------------------------------------Plot-------------------------------------------
prompt=sprintf('Please insert the desired atmittance as Ymn for m,n = 1,2,3 (e.g., Y11, Y12,...): ');
str1=input(prompt,'s');

str2={'Y11','Y12','Y13','Y14';'Y21','Y22','Y23', 'Y24';'Y31','Y32','Y33', 'Y34';'Y41','Y42','Y43', 'Y44'};
Ymn = strcmp(str1,str2);

if Ymn(1,1) == 1
    original=Y11_original;
    reduced=Y11_reduced;
elseif Ymn(1,2) == 1 || Ymn(2,1) == 1
    original=Y21_original;
    reduced=Y21_reduced;
elseif Ymn(1,3) == 1 || Ymn(3,1) == 1
    original=Y31_original;
    reduced=Y31_reduced;
elseif Ymn(1,4) == 1 || Ymn(4,1) == 1
    original=Y41_original;
    reduced=Y41_reduced;
elseif Ymn(2,2) == 1
    original=Y22_original;
    reduced=Y22_reduced;
elseif Ymn(2,3) == 1 || Ymn(3,2) == 1
    original=Y32_original;
    reduced=Y32_reduced;
elseif Ymn(2,4) == 1 || Ymn(4,2) == 1
    original=Y42_original;
    reduced=Y42_reduced;    
elseif Ymn(3,3) == 1
    original=Y33_original;
    reduced=Y33_reduced;    
elseif Ymn(3,4) == 1 || Ymn(4,3) == 1
    original=Y43_original;
    reduced=Y43_reduced; 
elseif Ymn(4,4) == 1
    original=Y44_original;
    reduced=Y44_reduced;
end

figure
semilogy(f,abs(original),'-.r','LineWidth',3)
hold on
semilogy(f,abs(reduced),'--g','LineWidth',3)
legend('MNA_{original}','MNA_{reduced}')
xlabel('f (THz)') 
ylabel(sprintf('| Y_{%s} |',str1(2:3))) 
title(sprintf('Y_{%s}',str1(2:3))) 

RMSE = sqrt(mean((abs(original)-abs(reduced)).^2))
title(sprintf('RMSE Reduced (n = %d) vs Original (N=50) = %d',n,RMSE))

grid on
hold on
clear Ymn str1 str2 prompt original reduced A B0 B1 Bm Bmr Cm Cmr Coef C D f fmax fmin freq fstep Gamma Gammar Gm Gmr H i I_orig11 I_orig21 I_orig31 I_orig22 I_orig32 I_orig33 I_red11 I_red21 I_red31 I_red22 I_red32 I_red33 j Jm k K K_inv L Lv Lvr M nt Num_branch p p0 q Q Q0 P0 QP QP0 QP_hat QPm Ra Rb s s0 V_orig V_reduced