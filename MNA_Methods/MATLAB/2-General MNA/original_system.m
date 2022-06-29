%---------------------------------------------------------------------------------------------------------------------------------------%

% original_system.m
% 
% Author: Rasul Choupanzadeh
% Date: 05/12/2022
% 
% 
% This code is based on the concepts from [1-2].
% 
% [1] Yangfeng Su, Jian Wang, Xuan Zeng, Zhaojun Bai, C. Chiang, and D. Zhou.
%     SAPOR: second-order arnoldi method for passive order reduction of RCS circuits.
%     In IEEE/ACM International Conference on Computer Aided Design (ICCAD), pages 74–79, 2004.
% 
% [2] Bang Liu, Xuan Zeng, Yangfeng Su, Jun Tao, Zhaojun Bai, C. Chiang, and Dian Zhou.
%     Block SAPOR: block second-order arnoldi method for passive order reduction of multi-input multi-output 
%     RCS interconnect circuits. In Proceedings of Asia and South Pacific Design Automation Conference (ASP-DAC), pages 244–249, 2005.

%---------------------------------------------------------------------------------------------------------------------------------------%



V_orig=zeros(q,length(f));
nt=1;
for freq=fmin:fstep:fmax
    s=1i*2*pi*freq*Coef;
    V_orig(:,nt)=transpose(Lv)*((s.*Cm+Gm+Gamma./s)\Bm)*Jm;
    nt=nt+1;
end 

clear nt freq s 