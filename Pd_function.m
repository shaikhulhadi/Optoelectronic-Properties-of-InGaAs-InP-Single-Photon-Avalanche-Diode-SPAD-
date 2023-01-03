function F=Pd_function(pd)
 global Pa;
%% list of given values 
tou = 2e-6;          %second ----- gate pulse width
tou_d=200e-9;        %characteristics detraping time constant
N_tr = 1e8;          % average number of trapped carrier  after a current pulse
N_t0 = N_tr/99;      % trapped carrier per pulse----- 1% of total carrier----- there is an equation
dekta_T=1/(100e3);   %time between successive gate pulse .....sum of tou and holdoff time
I_DM = 0.1e-12;     % primary dark current......1p /10p
M_o=20;             %DC gain....30/40
GB=30e9;            %gain bandwidth product
c=1;                %ration of trapped carries to the total carrier per avalence pulse
Mg=1e8;             %gaiger mode gain
frequency=100e3;    %frequency
q_charge=1.6e-19;   % charge
t_tr_star = M_o/(2*pi*GB);
                    %effective transit time
                    
%% copied from main code page
%% equation 
F=1-exp(-Pa*(I_DM*tou/q_charge+I_DM*M_o.^2/(2*pi*q_charge*GB)+...
    pd*N_tr*((exp(tou/tou_d)-1)/(exp(dekta_T/tou_d)-1)+...
    (exp(t_tr_star/tou_d)-1)/(exp(dekta_T/tou_d)-1))))-pd;
end
