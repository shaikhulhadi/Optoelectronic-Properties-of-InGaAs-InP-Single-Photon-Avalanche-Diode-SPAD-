%{
This Matlab code simulates darc count probability vs single
photon detection efficiency graph along with experimental value for 
-50,-30,0 and 20 degree celsius
%}
close all;
clear all;
clc;

%% list of parameter values
frequency=1e9;              %Gate frequency
delta_T=1/(frequency);      %time between successive gate pulse 
tou = 400e-12;              %in second ----- gate pulse width
tou_d=(delta_T-tou)/50 ;    %hold of time of the gate
electron_charge=1.6e-19;    % electron charge
N_0=0.10;                   %photon per pulse-----(+- 5%)
P_ph=1-exp(-N_0);           %probability of pulse containing photon
N_tr =1e8;                  % trapped carrier  per pulse
N_tr0 = 1e2;                % trapped carrier per pulse----- 1% of total carrier
I_DM = [.00045e-12 .001e-12 .001e-12    .1e-13];     % primary dark current......calibrated
M_o=[12 5 20 12];            %assumed
GB=30e9;                    %gain bandwidth product
c=0.01;                     %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;                     %gaiger mode gain

tolarence_pd=1e-10;         %tolarance values for itterative algorithom 
tolarence_pon=1e-10;

 %%experimental values reading from csv file
Pd_blue_50c=csvread('blue_50c.csv');  %% -50 degree celsius
Pd_b=Pd_blue_50c(:,2)';
SPDE_b=Pd_blue_50c(:,1)';

Pd_olive_30c=csvread('olive_30c.csv'); %%-30 degree celsius
Pd_o=Pd_olive_30c(:,2)';
SPDE_o=Pd_olive_30c(:,1)';

Pd_green_0c=csvread('green_0c.csv');    %%zero degree celsius
Pd_g=Pd_green_0c(:,2)';
SPDE_g=Pd_green_0c(:,1)';

Pd_red_20c=csvread('red_20c.csv');      %%20 degree celsius
Pd_r=Pd_red_20c(:,2)';
SPDE_r=Pd_red_20c(:,1)';


len=[length(Pd_b)   length(Pd_o)    length(Pd_g)    length(Pd_r)];
        %%%% lenght of vectors to avoid matrix mismatch error

for j=1:4
    t_tr_star = M_o(j)/(2*pi*GB);
    Pava=linspace(0.01,.8,len(j));  %% avalanche probability
    for i=1:1:length(Pava)
        Pa=Pava(i);
        delta_pd=1;
        pd_val=0.01;     
        delta_pon=1;
        pon_val=0.02;
        pd_temp=0;
        pon_temp=0;
        iteration(j,i) =0;  %% counting iteration value for debugging 
        QE=0.550;           %%greater QE yields better range of SPDE
        while (tolarence_pd < delta_pd) && ( tolarence_pon  <  delta_pon )
            if tolarence_pd<delta_pd
                pd_temp=1-exp(-Pa*(I_DM(j)*tou/electron_charge + ...
                    I_DM(j).*(M_o(j) ).^2/(2*pi*electron_charge*GB)+...
                    pd_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                    (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
                delta_pd=pd_val-pd_temp;
                pd_val=pd_temp;
            end
            
            if tolarence_pd < delta_pd
                pon_temp = 1-exp(-Pa*(I_DM(j)*tou/electron_charge + ...
                    I_DM(j)*M_o(j).^2/(2*pi*electron_charge*GB)+...
                    pon_val*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                    (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE*N_0));
                delta_pon = pon_val - pon_temp;
                pon_val=pd_temp;
            end
            SPDE_temp=(pon_temp - pd_temp)/P_ph;
            QE=SPDE_temp/Pa;
            iteration(j,i) =iteration(j,i) +1;
        end
        
        SPDE(j,i)=100*SPDE_temp;
        Pd(j,i)=pd_temp;
        Pon(j,i)=pon_temp;
    end
    
end
%%plot simulated results 
semilogy(SPDE(2,:),Pd(2,:),'m',SPDE(3,:),Pd(3,:),'g',SPDE(4,:),Pd(4,:),'r');
xlabel('Single-Photon Detection Efficiency %');
ylabel('Dark Count Probability,Pd'); 
legend('-30 degree','0 degree','20 degree','Location','southeast');
hold on
%% plot experimental data along with simulated data
semilogy(SPDE_o,Pd_o,'m--',SPDE_g,Pd_g,'g--',SPDE_r,Pd_r,'r--')