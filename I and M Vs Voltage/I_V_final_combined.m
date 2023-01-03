%{
This Matlab code simulates photocurrent vs applied voltage
graph along with experimental value for 
-50,-30,0 and 20 degree celsius
%}
close all
clear all
clc

Plank_constant = 6.63e-34;
Light_speed=3e8;
Electron_charge=1.6e-19;
Lambda = 1550e-9;
Input_Power = 196e-9;
n=5;                                    %% calibrated for better curve
Temp=[-50 -30 0 20];                    %% Temparatures
Temp_not=-50;
Vbro=60.1;                              %%breadown voltage at temperature Temp_not celsius 
V_punch=36 ;                            %%from the extracted values       
V_applied=30:.1:70;                     %%applied voltages
R = 0.3*(Electron_charge*Lambda)/(Plank_constant*Light_speed);   % responcivity 
Ipho = R*Input_Power;                   % unity gain current
Idm=1e-7*ones(1,length(V_applied)) ;    %%dark current 
index_punch=find(V_applied==V_punch); 
M=zeros(1,length(V_applied));  
M(index_punch)=1; 
%%%%%% photocurrent for different applied voltage

Vbr=Vbro+0.1*(20-Temp_not);     %% for 20 degree temperature
index_break=find(V_applied==Vbr);
T=V_applied/Vbr;
for i=index_punch+1:1:index_break 
 M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)))+(((T(i))^(10*n)))+...
            (((T(i))^(11*n)))+(((T(i))^(12*n)))+(((T(i))^(13*n)))+(((T(i))^(14*n)))+(((T(i))^(15*n)))...
            +(((T(i))^(16*n)));
end
for i=index_break+1:1:length(V_applied)
    M(i)=1/0;
end
Id=M.*Idm+.02*Idm;
Iout=M.*Ipho+Id;

 red_20c=csvread('red.csv');  %%%% data extraction and plotting for -20 degree celsius
Iph_exp_r(1,:)=red_20c(:,2)';
V_exp_r(1,:)=red_20c(:,1)';
simu_data(:,7)=V_applied;
simu_data(:,8)=Iout;
plot(V_applied,Iout,'r');


Vbr=Vbro+0.1*(0-Temp_not);     %% for 0 degree temperature
index_break=find(V_applied==Vbr);
T=V_applied/Vbr;
for i=index_punch+1:1:index_break
 M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)))+(((T(i))^(10*n)))+...
            (((T(i))^(11*n)))+(((T(i))^(12*n)))+(((T(i))^(13*n)))+(((T(i))^(14*n)))+(((T(i))^(15*n)))...
            +(((T(i))^(16*n)));
end
for i=index_break+1:1:length(V_applied)
    M(i)=1/0;
end
Id=M.*Idm+.02*Idm;
Iout=M.*Ipho+Id;
green_0c=csvread('green.csv');  %% experimental data for zero degree celsius
Iph_exp_g=green_0c(:,2)';
V_exp_g=green_0c(:,1)';
hold on
simu_data(:,5)=V_applied;
simu_data(:,6)=Iout;
plot(V_applied,Iout,'g');
legend( '0 degree' )

Vbr=Vbro+0.1*(-30-Temp_not);     %% for -30 degree temperature
index_break=find(V_applied==Vbr);
T=V_applied/Vbr;
for i=index_punch+1:1:index_break
 M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)))+(((T(i))^(10*n)))+...
            (((T(i))^(11*n)))+(((T(i))^(12*n)))+(((T(i))^(13*n)))+(((T(i))^(14*n)))+(((T(i))^(15*n)))...
            +(((T(i))^(16*n)));
end
for i=index_break+1:1:length(V_applied)
    M(i)=1/0;
end
Id=M.*Idm+.02*Idm;
Iout=M.*Ipho+Id;
olive_30c=csvread('olive.csv');
Iph_exp_o=olive_30c(:,2)';
V_exp_o=olive_30c(:,1)';
hold on;
simu_data(:,3)=V_applied;
simu_data(:,4)=Iout;
plot(V_applied,Iout,'m');
legend('-30 degree');


Vbr=Vbro+0.1*(-50-Temp_not);     %% for -50 degree temperature
index_break=find(V_applied==Vbr);
T=V_applied/Vbr;
for i=index_punch+1:1:index_break
 M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)))+(((T(i))^(10*n)))+...
            (((T(i))^(11*n)))+(((T(i))^(12*n)))+(((T(i))^(13*n)))+(((T(i))^(14*n)))+(((T(i))^(15*n)))...
            +(((T(i))^(16*n)));
end
for i=index_break+1:1:length(V_applied)
    M(i)=1/0;
end
Id=M.*Idm+.02*Idm;
Iout=M.*Ipho+Id;
blue_50=csvread('blue.csv');
Iph_exp_b=blue_50(:,2)';
V_exp_b=blue_50(:,1)';
simu_data(:,1)=V_applied;
simu_data(:,2)=Iout;
hold on
plot(V_applied,Iout,'b');
hold on
plot(V_exp_r,Iph_exp_r,'r--',V_exp_g,Iph_exp_g,'g--',V_exp_o,Iph_exp_o,'m--',V_exp_b,Iph_exp_b,'b--')
legend('20 degree','0 degree','-30 degree','-50 degree','Location','northwest');
xlabel('Voltage,V');
ylabel('Photo current, Ip');
csvwrite('simu_data.csv',simu_data);  %% value is saved to plot graph in origin


 