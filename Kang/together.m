close all
clear all
clc


%% list of given values
frequency=100e3;    %frequency`````````````````````````````````````````````````````````````````````````(C)
delta_T=1/(frequency);  %time between successive gate pulse .....sum of tou and holdoff time``````````(C)
tou = 2e-9;     % delta_t/5000    %second ----- gate pulse width````````````````````````````````````````````````(C)
tou_d=200e-9;; %hold of time of the gate............maximum detrapping time ```````````````(C)....calculated
%tou_d,smaller the better
q_charge=1.6e-19;   % charge
%effective transit time
N_0=0.30;            %photon per pulse`````````````````````````````````````````````````(C)+- 5%
P_ph=1-exp(-N_0);     %probability of pulse containing photon

%---------------------------------------%
N_tr =1e8;          % trapped carrier  per pulse
N_tr0 = 1e6;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
I_DM = .1e-12;     % primary dark current......1p /10p....kang
M_o=10;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=0.01;                %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;             %gaiger mode gain
t_tr_star = M_o/(2*pi*GB);
%---------------------------------------%
%QE = 0.8;           %quantum efficiency...eta in commander .69


%P_ph=1-exp(-.1)    %probability of pulse containing photon





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations%%%%%%%%%%%%%%%%
tolarence_pd=1e-10;
tolarence_pon=1e-10;

%%%% Readin the values from cvs file
Pd_01pA=csvread('experimental_Idm_0.1pA.csv');
Pd_expe=Pd_01pA(:,2)';
SPDE_expe=Pd_01pA(:,1)';
Pd_01pA_sim=csvread('simulation_Idm_0.1pA.csv');
Pd_sim=Pd_01pA_sim(:,2)';
SPDE_sim=Pd_01pA_sim(:,1)';

Pava=linspace(0.01,1,length(Pd_expe));    %avalence probability assumed 
%%Avalanche probability range is 
%very important parameter
%Avalanche probabilit  yiels better range extention for SPDE
%%% linspace is used so that number of Pd and number of Pd_experiment data
%%% match to observe the difference in curve
%SPDE=1;
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
    QE=0.60;%%QE greater yields better range of SPDE .7 is very good
    while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
        if tolarence_pd<delta_pd
            pd_temp=1-exp(-Pa*(I_DM*tou/q_charge + I_DM.*(M_o ).^2/(2*pi*q_charge*GB)+...
                pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
            delta_pd=pd_val-pd_temp;
            delta_pd_sum=delta_pd_sum+delta_pd;  %
            pd_val=pd_temp;
        end
        
        if tolarence_pon < delta_pon
            pon_temp = 1-exp(-Pa*(I_DM*tou/q_charge + I_DM*M_o.^2/(2*pi*q_charge*GB)+...
                pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
            delta_pon = pon_val - pon_temp;
            delta_pon_sum=delta_pon_sum+delta_pon;%
            pon_val=pon_temp;
        end
        SPDE_temp=(pon_temp - pd_temp)/P_ph;
        QE=SPDE_temp/Pa;
        iteration(i) =iteration(i)+1;
    end
    QE_all(i)=QE;
    SPDE(i)=SPDE_temp;
    SPDE(i)=100*SPDE(i);
    Pd(i)=pd_temp;
    Pon(i)=pon_temp;
end
figure(1)
semilogy(SPDE,Pd,'k',SPDE_expe,Pd_expe,'r',SPDE_sim,Pd_sim,'g');
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability Per Pulse'); 
legend('Our Simulation ','Kang Experimental','Kang Simulation');






close all
clear all
clc

%% list of given values

%% list of given values
frequency=100e3;    %frequency`````````````````````````````````````````````````````````````````````````(C)
delta_T=1/(frequency);  %time between successive gate pulse .....sum of tou and holdoff time``````````(C)
tou = 2e-9;     % delta_t/5000    %second ----- gate pulse width````````````````````````````````````````````````(C)
tou_d=200e-9;; %hold of time of the gate............maximum detrapping time ```````````````(C)....calculated
%tou_d,smaller the better
q_charge=1.6e-19;   % charge
%effective transit time
N_0=0.30;            %photon per pulse`````````````````````````````````````````````````(C)+- 5%
P_ph=1-exp(-N_0);     %probability of pulse containing photon

%---------------------------------------%
N_tr =1e8;          % trapped carrier  per pulse
N_tr0 = 1e6;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
I_DM = 1e-12;     % primary dark current......1p /10p....kang
M_o=20;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=0.01;                %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;             %gaiger mode gain
t_tr_star = M_o/(2*pi*GB);
%---------------------------------------%
%QE = 0.8;           %quantum efficiency...eta in commander .69


%P_ph=1-exp(-.1)    %probability of pulse containing photon





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations%%%%%%%%%%%%%%%%
tolarence_pd=1e-10;
tolarence_pon=1e-10;

%%%% Readin the values from cvs file
Pd_1pA=csvread('Experimental_Idm_1pA.csv');
Pd_expe=Pd_1pA(:,2)';
SPDE_expe=Pd_1pA(:,1)';

Pd_1pA=csvread('simulation_Idm_1pA.csv');
Pd_sim=Pd_1pA(:,2)';
SPDE_sim=Pd_1pA(:,1)';
Pava=linspace(0.01,1,length(Pd_expe));    %avalence probability assumed 
%%Avalanche probability range is 
%very important parameter
%Avalanche probabilit  yiels better range extention for SPDE
%%% linspace is used so that number of Pd and number of Pd_experiment data
%%% match to observe the difference in curve
%SPDE=1;
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
    QE=0.6;%%QE greater yields better range of SPDE .7 is very good
    while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
        if tolarence_pd<delta_pd
            pd_temp=1-exp(-Pa*(I_DM*tou/q_charge + I_DM.*(M_o ).^2/(2*pi*q_charge*GB)+...
                pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
            delta_pd=pd_val-pd_temp;
            delta_pd_sum=delta_pd_sum+delta_pd;  %
            pd_val=pd_temp;
        end
        
        if tolarence_pon < delta_pon
            pon_temp = 1-exp(-Pa*(I_DM*tou/q_charge + I_DM*(M_o).^2/(2*pi*q_charge*GB)+...
                pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
            delta_pon = pon_val - pon_temp;
            delta_pon_sum=delta_pon_sum+delta_pon;%
            pon_val=pon_temp;
        end
        SPDE_temp=(pon_temp - pd_temp)/P_ph;
        QE=SPDE_temp/Pa;
        iteration(i) =iteration(i)+1;
    end
    SPDE(i)=SPDE_temp;
    SPDE(i)=100*SPDE(i);
    Pd(i)=pd_temp;
    Pon(i)=pon_temp;
end
figure(1)
semilogy(SPDE,Pd,'k',SPDE_expe,Pd_expe,'r',SPDE_sim,Pd_sim,'g');
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability Per Pulse'); 
legend('Our Simulation ','Kang Experimental','Kang Simulation');

%xlim([0 60]);
%  figure(2)
%  subplot(2,1,1)%
% semilogy(SPDE,Pd,'k*');
%  xlabel('Single-Photon Detection Efficiency %');
%  ylabel('Dark Count Probability Per Pulse'); 
%  
%  
%  subplot(2,1,2)
%  semilogy(SPDE_expe,Pd_expe,'r*');
% xlabel('Single-Photon Detection Efficiency %');
%  ylabel('Dark Count Probability Per Pulse'); 





close all
clear all
clc

%% list of given values

%% list of given values
frequency=100e3;    %frequency`````````````````````````````````````````````````````````````````````````(C)
delta_T=1/(frequency);  %time between successive gate pulse .....sum of tou and holdoff time``````````(C)
tou = 2e-9;     % delta_t/5000    %second ----- gate pulse width````````````````````````````````````````````````(C)
tou_d=200e-9;; %hold of time of the gate............maximum detrapping time ```````````````(C)....calculated
%tou_d,smaller the better
q_charge=1.6e-19;   % charge
%effective transit time
N_0=0.30;            %photon per pulse`````````````````````````````````````````````````(C)+- 5%
P_ph=1-exp(-N_0);     %probability of pulse containing photon

%---------------------------------------%
N_tr =1e8;          % trapped carrier  per pulse
N_tr0 = 1e6;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
I_DM = 10e-12;     % primary dark current......1p /10p....kang
M_o=10;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=0.01;                %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;             %gaiger mode gain
t_tr_star = M_o/(2*pi*GB);
%---------------------------------------%
%QE = 0.8;           %quantum efficiency...eta in commander .69


%P_ph=1-exp(-.1)    %probability of pulse containing photon





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations%%%%%%%%%%%%%%%%
tolarence_pd=1e-10;
tolarence_pon=1e-10;

%%%% Readin the values from cvs file
Pd_10pA=csvread('Theoratical_Idm_10pA.csv');
Pd_theo=Pd_10pA(:,2)';
SPDE_theo=Pd_10pA(:,1)';





Pava=linspace(0.01,1,length(Pd_theo));    %avalence probability assumed 
%%Avalanche probability range is 
%very important parameter
%Avalanche probabilit  yiels better range extention for SPDE
%%% linspace is used so that number of Pd and number of Pd_experiment data
%%% match to observe the difference in curve
%SPDE=1;
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
    QE=0.60;%%QE greater yields better range of SPDE .7 is very good
    while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
        if tolarence_pd<delta_pd
            pd_temp=1-exp(-Pa*(I_DM*tou/q_charge + I_DM.*(M_o ).^2/(2*pi*q_charge*GB)+...
                pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
            delta_pd=pd_val-pd_temp;
            delta_pd_sum=delta_pd_sum+delta_pd;  %
            pd_val=pd_temp;
        end
        
        if tolarence_pon < delta_pon
            pon_temp = 1-exp(-Pa*(I_DM*tou/q_charge + I_DM*M_o.^2/(2*pi*q_charge*GB)+...
                pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
            delta_pon = pon_val - pon_temp;
            delta_pon_sum=delta_pon_sum+delta_pon;%
            pon_val=pon_temp;
        end
        SPDE_temp=(pon_temp - pd_temp)/P_ph;
        QE=SPDE_temp/Pa;
        iteration(i) =iteration(i)+1;
    end
    SPDE(i)=SPDE_temp;
    SPDE(i)=100*SPDE(i);
    Pd(i)=pd_temp;
    Pon(i)=pon_temp;
end
figure(1)
semilogy(SPDE,Pd,'k',SPDE_theo,Pd_theo,'r');
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability Per Pulse'); 
legend('Our Simulation ','Kang Simulation or theoretical ');

%  figure(2)
%  subplot(2,1,1)%
% semilogy(SPDE,Pd,'k*');
%  xlabel('Single-Photon Detection Efficiency %');
%  ylabel('Dark Count Probability Per Pulse'); 
%  
%  
%  subplot(2,1,2)
%  semilogy(SPDE_expe,Pd_expe,'r*');
% xlabel('Single-Photon Detection Efficiency %');
%  ylabel('Dark Count Probability Per Pulse'); 