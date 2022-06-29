%---------------------------------------------------------%

% Load_netlist.m
% 
% Author: Rasul Choupanzadeh
% Date: 05/12/2022

%---------------------------------------------------------%




function [Ra,Rb,L,C,Num_branch]=Load_netlist(sbck)
global netlist ;
file='netlist_passive.sp';
fileID = fopen(file);
netlist = textscan(fileID,'%s %s %s %s %s %s %s %s %s');
fclose(fileID);
clear fileID file ans;
%-------------------------------Split the sub-circuits------------------------------------
subckt_identifier=[];
nt=1;
for i=1:length(netlist{1})
    s = netlist{1}{i};
    st= '.ends';
    tf = strcmp(s,st);                          % Compares the strings with '.ends'
    if tf==1
        subckt_identifier(nt)=i+1;              % Stores the row numbers for the begining of each sub_circuit; which is occures after each '.ends'
        nt=nt+1;
    end  
end

% Extracting parameters of sub-circuit to calculate MNA matrices
Num_sbckt=length(subckt_identifier)-1; 
Num_branch=0;
for i=subckt_identifier(sbck):subckt_identifier(sbck+1)-1
    s = netlist{2}{i};
    st= 'Branch';
    tf = strcmp(s,st);                      
    if tf==1
        Num_branch=Num_branch+1;              
    end  
end
Ra=zeros(1,Num_branch);
L=zeros(1,Num_branch);
Rb=zeros(1,Num_branch);
C=zeros(1,Num_branch);
nt1=1;
nt2=1;
nt3=1;
nt4=1;
for i=subckt_identifier(sbck):subckt_identifier(sbck+1)-1
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
end


