%---------------------------------------------------------------------------------------------------------------------------------------%

% Block_SAPOR.m
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


% Prepating input matrices
D=2*s0*Cm+Gm;
K=(s0^2)*Cm+s0*Gm+Gamma;
B0=s0*Bm;
B1=Bm;
K_inv=K\eye(length(K));
A=[-K_inv*D K_inv;-Cm zeros(N,N)];
Q0=K_inv*B0;    
P0=B1;

% Block SAPOR main algorithm (inputs: [A],[Q0],[P0],n,p     output: projection matrix [Q])
QP0=[Q0;P0];
QP_hat=[];
k=n/p;
[QPm]=SOrth(QP0,p,N); 
QP=QPm;
for i=1:k-1
    QP_hat(:,i*p+1:(i+1)*p)=A*QP(:,(i-1)*p+1:i*p);
    for j=1:i
        H=transpose(QP(1:N,(j-1)*p+1:j*p))*QP_hat(1:N,i*p+1:(i+1)*p);
        QP_hat(:,i*p+1:(i+1)*p)=QP_hat(:,i*p+1:(i+1)*p)-QP(:,(j-1)*p+1:j*p)*H;
    end
    [QPm]=SOrth(QP_hat(:,i*p+1:(i+1)*p),p,N);
    QP(:,i*p+1:(i+1)*p)=QPm;
end
Q=QP(1:N,:);

% Calcualting reduced MNA matrices projected by [Q]
Cmr=transpose(Q)*Cm*Q;
Gmr=transpose(Q)*Gm*Q;
Gammar=transpose(Q)*Gamma*Q;
Bmr=transpose(Q)*Bm;
Lvr=transpose(Q)*Lv;

% Calculating frequency response
tic
V_reduced=zeros(q,length(f));
nt=1;
for freq=fmin:fstep:fmax 
    s=1i*2*pi*freq*Coef;
    V_reduced(:,nt)=transpose(Lvr)*((s.*Cmr+Gmr+Gammar./s)\Bmr)*Jm;
    nt=nt+1;
end
toc
%clear D K B0 B1 K_inv A k nt freq s Q0 P0 QP_hat QP QP0 QPm s0 i j H fmax fmin fstep M 