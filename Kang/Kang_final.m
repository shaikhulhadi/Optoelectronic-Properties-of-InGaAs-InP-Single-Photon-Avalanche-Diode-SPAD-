close all
clear all
clc

%% list of given values
frequency=100e3;                        %frequency
delta_T=1/(frequency);          %time between successive gate pulse
tou = 2e-9;                     %second ----- gate pulse width
tou_d=200e-9;                   %hold of time of the gate
q_charge=1.6e-19;               % charge
N_0=0.30;                       %photon per pulse(+- 5%)
P_ph=1-exp(-N_0);               %probability of pulse containing photon
N_tr =1e8;                      % trapped carrier  per pulse
N_tr0 = 1e6;                    % trapped carrier per pulse
I_DM = [.1e-12  1e-12   10e-12]; % primary dark current......1p /10p....kang
M_o=10;                         %DC gain
GB=30e9;                        %gain bandwidth product
c=0.01;                         %ratio of trapped carries to the total carrier per avalanche pulse
Mg=1e8;                         %gaiger mode gain
t_tr_star = M_o/(2*pi*GB);
tolarence_pd=1e-10;
tolarence_pon=1e-10;

%%%% Readin the values from cvs file
Pd_01pA=csvread('experimental_Idm_0.1pA.csv');
Pd_ex_01=Pd_01pA(:,2)';
SPDE_ex_01=Pd_01pA(:,1)';
Pd_01pA_sim=csvread('simulation_Idm_0.1pA.csv');
Pd_sim_01=Pd_01pA_sim(:,2)';
SPDE_sim_01=Pd_01pA_sim(:,1)';

Pd_1pA=csvread('Experimental_Idm_1pA.csv');
Pd_ex_1=Pd_1pA(:,2)';
SPDE_ex_1=Pd_1pA(:,1)';
Pd_1pA=csvread('simulation_Idm_1pA.csv');
Pd_sim_1=Pd_1pA(:,2)';
SPDE_sim_1=Pd_1pA(:,1)';

Pd_10pA=csvread('Theoratical_Idm_10pA.csv');
Pd_th_10=Pd_10pA(:,2)';
SPDE_th_10=Pd_10pA(:,1)';

len=[length(Pd_ex_01)   length(Pd_ex_1)    length(Pd_th_10)];

for j=1:3
    Pava=linspace(0.01,.7,len(j));
    for i=1:1:length(Pava)
        Pa=Pava(i);
        delta_pd=1;
        pd_val=0.01;
        
        delta_pon=1;
        pon_val=0.02;
        pd_temp=0;
        pon_temp=0;
        iteration(i) =0;
        QE=0.60;
        while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
            if tolarence_pd<delta_pd
                pd_temp=1-exp(-Pa*(I_DM(j)*tou/q_charge + I_DM(j).*(M_o ).^2/(2*pi*q_charge*GB)+...
                    pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                    (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
                delta_pd=pd_val-pd_temp;
                pd_val=pd_temp;
            end
            
            if tolarence_pon < delta_pon
                pon_temp = 1-exp(-Pa*(I_DM(j)*tou/q_charge + I_DM(j)*M_o.^2/(2*pi*q_charge*GB)+...
                    pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                    (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
                delta_pon = pon_val - pon_temp;
                pon_val=pon_temp;
            end
            SPDE_temp=(pon_temp - pd_temp)/P_ph;
            QE=SPDE_temp/Pa;
            iteration(i) =iteration(i)+1;
        end
     
        SPDE(j,i)=100*SPDE_temp;
        Pd(j,i)=pd_temp;
        Pon(j,i)=pon_temp;
    end
end
semilogy(SPDE(1,:),Pd(1,:),'b',SPDE(2,:),Pd(2,:),'g',SPDE(3,:),Pd(3,:),'r');
            %% simulation data plot
hold on
semilogy(SPDE_ex_01,Pd_ex_01,'b--',SPDE_ex_1,Pd_ex_1,'g--',SPDE_sim_01,Pd_sim_01,'b+'...
    ,SPDE_sim_1,Pd_sim_1,'g+',SPDE_th_10,Pd_th_10,'r+');
            %%experimental data plot
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability Per Pulse'); 
legend('0.1pA ','1 pA','10 pA','Location','southeast');
