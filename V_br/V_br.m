close all
clear all
clc

V_bro = 60;
T_o=-50;
T_app=-60:1:25;
gamma=0.1;
V_bre=V_bro +gamma.*(T_app-T_o);
plot(V_bre,T_app)
for i=1:length(T_app)
data(i,1)=T_app(i);
data(i,2)=V_bre(i);
end 

xlswrite('ank.xlsx',data);