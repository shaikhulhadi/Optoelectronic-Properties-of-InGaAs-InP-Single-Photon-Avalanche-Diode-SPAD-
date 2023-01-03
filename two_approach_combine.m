close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% list of given values

%% list of given values
frequency=1e9;    %frequency`````````````````````````````````````````````````````````````````````````(C)
delta_T=1/(frequency);   %time between successive gate pulse .....sum of tou and holdoff time``````````(C)
tou = 360e-12;          %second ----- gate pulse width````````````````````````````````````````````````(C)
q_charge=1.6e-19;   % charge
%effective transit time
N_0=0.1;            %photon per pulse`````````````````````````````````````````````````(C)+- 5%


%---------------------------------------%
N_tr = 1e8;          % trapped carrier  per pulse
N_tr0 = 1e6;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
I_DM = 1e-13;     % primary dark current......1p /10p....kang
M_o=20;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=0.01;                %ration of trapped carries to the total carrier per avalanche pulse
Mg=1e8;             %gaiger mode gain
tou_d=(delta_T-tou)/200;        %hold of time of the gate............maximum detrapping time ```````````````....calculated
t_tr_star = M_o/(2*pi*GB);
%---------------------------------------%

P_ph=1-exp(-N_0);     %probability of pulse containing photon

%%%% Readin the values from cvs file
Pd_red_20c=csvread('red_20c.csv');
Pd_curve=Pd_red_20c(:,2)';
SPDE_curve=100*Pd_red_20c(:,1)'/100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pava=linspace(0.0003,0.09,length(Pd_curve));    %avalence probability assumed
tolarence_pd=1e-10;
tolarence_pon1=1e-10;


%%%% equation to calculate Pa from Pd through approach 1


for i=1:1:length(Pava)
    Pa1=Pava(i);
    delta_pd1=1;
    delta_pd_old1=0;
    delta_pd_sum1=0;
    pd_val1=0.01;
    
    delta_pon1=1;
    delta_pon_old2=0;
    delta_pon_sum1=0;
    pon_val1=0.02;
    pd_temp1=0;
    pon_temp1=0;
    iteration1(i) =0;
    QE1=0.8;
    while (tolarence_pd < delta_pd1) && ( tolarence_pon1  <  delta_pon1 )
        if tolarence_pd<delta_pd1
            pd_temp1=1-exp(-Pa1*(I_DM*tou/q_charge + I_DM.*(M_o ).^2/(2*pi*q_charge*GB)+...
                pd_val1*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1))));
            delta_pd1=pd_val1-pd_temp1;
            delta_pd_sum1=delta_pd_sum1+delta_pd1;  %
            pd_val1=pd_temp1;
        end
        
        if tolarence_pd < delta_pd1
            pon_temp1 = 1-exp(-Pa1*(I_DM*tou/q_charge + I_DM*M_o.^2/(2*pi*q_charge*GB)+...
                pon_val1*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
                (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE1*N_0));
            delta_pon1 = pon_val1 - pon_temp1;
            delta_pon_sum1=delta_pon_sum1+delta_pon1;  %
            pon_val1=pd_temp1;
        end
        SPDE_temp1=(pon_temp1 - pd_temp1)/P_ph;
        QE1=SPDE_temp1/Pa1
        iteration1(i) =iteration1(i) +1;
    end
    SPDE(i)=100*SPDE_temp1;
    Pd1(i)=pd_temp1;
    Pon1(i)=pon_temp1;
end


%%%%%%%%%%%%% approach 2
Pd2= Pd_curve;
Pa = -log(1-Pd2)./(I_DM*tou/q_charge + I_DM.*(M_o ).^2./(2*pi*q_charge*GB)+...
    Pd2.*N_tr*((exp(tou/tou_d)-1)./(exp(delta_T/tou_d)-1)+...
    (exp(t_tr_star/tou_d)-1)./(exp(delta_T/tou_d)-1)));

%%%%%%%%% Now calculate the value of Pon through iteration


delta_pon2=ones(1,length(Pa));
delta_pon_old2=zeros(1,length(Pa));
delta_pon_sum2 = zeros(1,length(Pa));
Pon_val2=ones(1,length(Pa));
Pon_temp2=ones(1,length(Pa));
tolarence_pon2=1e-10;
iteration2 =0;
QE2=ones(1,length(Pa))*0.8

%%%%%%% iteration process to solve Pon and SPDE
while tolarence_pon2 < max(delta_pon2)
    for i=1:length(Pa)
        Pon_temp2(i) = 1-exp(-Pa(i)*(I_DM*tou/q_charge + I_DM*M_o.^2/(2*pi*q_charge*GB)+...
            Pon_val2(i)*N_tr*((exp(tou/tou_d)-1)/(exp(delta_T/tou_d)-1)+...
            (exp(t_tr_star/tou_d)-1)/(exp(delta_T/tou_d)-1)) + QE2(i).*N_0));
    end
    delta_pon2 = Pon_val2 - Pon_temp2;
    delta_pon_sum2=delta_pon_sum2+delta_pon2;  %
    Pon_val2=Pon_temp2;
    SPDE_temp2=(Pon_temp2 - Pd2)./P_ph;
    QE2=SPDE_temp2./Pa;
    iteration2 =iteration2 +1;
end
SPDE2=100*SPDE_temp2;
%figure(1)
%semilogy(Pd1,SPDE,'g',Pd2,SPDE2,'b',Pd_curve,SPDE_curve,'r')

figure(2)
subplot (2,2,1) 
semilogy(SPDE,Pd1,'g')
xlabel('SPDE');
ylabel('Dark Count Probability');
title('Approach One By Direct Formula and Iteration');
subplot(2,2,2) 
semilogy(SPDE2,Pd2,'b')
xlabel('SPDE');
ylabel('Dark Count Probability');
title('Approach Two,by back calculation Method');
subplot(2,2,3)
semilogy(SPDE_curve,Pd_curve,'r');
xlabel('SPDE');
ylabel('Dark Count Probability');
title('Graph from  L C Commander');

%semilogy(SPDE,Pd1,'g',SPDE2,Pd2,'b',SPDE_curve,Pd_curve,'r')



