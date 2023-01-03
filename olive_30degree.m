close all;
clear all;
clc;

%% list of given values

%% list of given values
frequency=1e9;    %frequency`````````````````````````````````````````````````````````````````````````(C)
delta_T=1/(frequency);  %time between successive gate pulse .....sum of tou and holdoff time``````````(C)
tou = 400e-12;     % delta_t/5000    %second ----- gate pulse width````````````````````````````````````````````````(C)
tou_d=(delta_T-tou)/200 ;  %hold of time of the gate............maximum detrapping time ```````````````(C)....calculated
q_charge=1.6e-19;   % charge
%effective transit time
N_0=0.10;            %photon per pulse`````````````````````````````````````````````````(C)+- 5%
P_ph=1-exp(-N_0);     %probability of pulse containing photon

%---------------------------------------%
N_tr =1e8;          % trapped carrier  per pulse
N_tr0 = 1e2;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
I_DM = .001e-12;     % primary dark current......1p /10p....kang
M_o=5;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=0.01;                %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;             %gaiger mode gain
t_tr_star = M_o/(2*pi*GB);
%---------------------------------------%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations%%%%%%%%%%%%%%%%
tolarence_pd=1e-10;
tolarence_pon=1e-10;

%%%% Readin the values from cvs file
Pd_red_20c=csvread('olive_30c.csv');
Pd_expe=Pd_red_20c(:,2)';
SPDE_expe=Pd_red_20c(:,1)';

Pava=linspace(0.01,.8,length(Pd_expe));    %avalence probability assumed 
%%Avalanche probability range is 
%very important parameter
%Avalanche probabilit  yiels better range extention for SPDE
%%% linspace is used so that number of Pd and number of Pd_experiment data
%%% match to observe the difference in curve

for i=1:1:length(Pava)
    Pa=Pava(i);
    delta_pd=1;
    delta_pd_old=0;
    delta_pd_sum=0;
    pd_val=0.01;
    
    delta_pon=1;
    delta_pon_old=0;
    delta_pon_sum=0;
    pon_val=0.02;
    pd_temp=0;
    pon_temp=0;
    iteration(i) =0;
    QE=0.550;%%QE greater yields better range of SPDE .7 is very good
    while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
        if tolarence_pd<delta_pd
            pd_temp=1-exp(-Pa*(I_DM*tou/q_charge + I_DM.*(M_o ).^2/(2*pi*q_charge*GB)+...
                pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
            delta_pd=pd_val-pd_temp;
            delta_pd_sum=delta_pd_sum+delta_pd;  %
            pd_val=pd_temp;
        end
        
        if tolarence_pd < delta_pd
            pon_temp = 1-exp(-Pa*(I_DM*tou/q_charge + I_DM*M_o.^2/(2*pi*q_charge*GB)+...
                pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
            delta_pon = pon_val - pon_temp;
            delta_pon_sum=delta_pon_sum+delta_pon;  %
            pon_val=pd_temp;
        end
        SPDE_temp=(pon_temp - pd_temp)/P_ph
        QE=SPDE_temp/Pa;
        iteration(i) =iteration(i) +1;
    end
    SPDE(i)=SPDE_temp;
    SPDE(i)=100*SPDE(i);
    Pd(i)=pd_temp;
    Pon(i)=pon_temp;
end
figure(1)
semilogy(SPDE,Pd,'k',SPDE_expe,Pd_expe,'r');
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability Per Pulse'); 

% figure(2)
% subplot(2,1,1)
% semilogy(SPDE,Pd,'k*');
% xlabel('Single-Photon Detection Efficiency %');
% ylabel('Dark Count Probability Per Pulse'); 
% 
% 
% subplot(2,1,2)
% semilogy(SPDE_expe,Pd_expe,'r*');