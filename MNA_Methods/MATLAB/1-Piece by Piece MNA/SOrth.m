%---------------------------------------------------------------------------------------------------------------------------------------%

% SOrth.m
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



function [QPm]=SOrth(QPm_hat,p,N)
for i=1:p
    qp(:,i)=QPm_hat(:,i);
    for j=1:i-1
        R(j,i)=transpose(qp(1:N,j))*qp(1:N,i);
        qp(:,i)=qp(:,i)-R(j,i)*qp(:,j);
    end
    R(i,i)=norm(qp(1:N,i));
    if R(i,i)==0
            break
    else
        qp(:,i)=R(i,i)\qp(:,i);
    end
end
QPm=qp;

