clear all
close all
clc
blue=1;
olive=2;
green=3;
red=4;
%%%%%%%%%%%%-------------------%%%%%%%%
Pd_blue_50c=csvread('blue_50c.csv',0,0,[0 0 34 1]);  %% number 1
Pd_b=Pd_blue_50c(:,2)';
SPDE_b=Pd_blue_50c(:,1)';

Pd_olive_30c=csvread('olive_30c.csv',0,0,[0 0 43 1]);  %% number 2
Pd_o=Pd_olive_30c(:,2)';
SPDE_o=Pd_olive_30c(:,1)';

Pd_green_0c=csvread('green_0c.csv',0,0,[0 0 34 1]);   %% number 3
Pd_g=Pd_green_0c(:,2)';
SPDE_g=Pd_green_0c(:,1)';

Pd_red_20c=csvread('red_20c.csv',0,0,[0 0 35 1]);  %% number 4
Pd_r=Pd_red_20c(:,2)';
SPDE_r=Pd_red_20c(:,1)';
%%%%%%%------------------%%%%%%%%%% simulated values from previous coding

SPDE_simu = csvread('SPDE_SIMULATED.csv');
Pd_simu=csvread('Pd_simulation.csv');

%%%%___________________%%%%%




%%%%%%%%---------calculation of  interpolation points due to  point mismatch --------%%%%%%

    %% ensuring no monotonicity which will cause error in interpolation algorithm 
    for i=1:length(SPDE_simu)
      if SPDE_simu(blue,i) ~= 0
          SPDE_b_simu(i)=SPDE_simu(blue,i);
          Pd_b_simu(i)=Pd_simu(blue,i);
      end
      if SPDE_simu(olive,i) ~= 0
          SPDE_o_simu(i)=SPDE_simu(olive,i);
          Pd_o_simu(i)=Pd_simu(olive,i);
      end
      if SPDE_simu(green,i) ~= 0
          SPDE_g_simu(i)=SPDE_simu(green,i);
          Pd_g_simu(i)=Pd_simu(green,i);
      end
      if SPDE_simu(red,i) ~= 0
          SPDE_r_simu(i)=SPDE_simu(red,i);
          Pd_r_simu(i)=Pd_simu(red,i);
      end
    end

    Pd_b_int= interp1(SPDE_b_simu,Pd_b_simu, SPDE_b); 
    Pd_o_int= interp1(SPDE_o_simu,Pd_o_simu, SPDE_o); 
    Pd_g_int= interp1(SPDE_g_simu,Pd_g_simu, SPDE_g); 
    Pd_r_int= interp1(SPDE_r_simu,Pd_r_simu, SPDE_r); 
    
  %%%% standard deviation of error percentage within 30 percent efficiency   
    stand_error(olive,1:length(Pd_o))=((Pd_o_int-Pd_o)./Pd_o*100);
    stand_error(blue,1:length(Pd_b))=((Pd_b_int-Pd_b)./Pd_b*100);
    stand_error(green,1:length(Pd_g))=((Pd_g_int-Pd_g)./Pd_g*100);
    stand_error(red,1:length(Pd_r))=((Pd_r_int-Pd_r)./Pd_r*100);
    
    
    
    %%%%%%%%%%%%%%%%
   %% SPDE=zeros(4,length(SPDE_o));
   
    SPDE_data(blue,1:length(SPDE_b))=SPDE_b;
    SPDE_data(olive,1:length(SPDE_o))=SPDE_o;
    SPDE_data(green,1:length(SPDE_g))=SPDE_g;
    SPDE_data(red,1:length(SPDE_r))=SPDE_r;
    
%     Error(blue,1:length(Pd_b))=abs((Pd_b_int-Pd_b))./Pd_b*100;
%     Error(olive,1:length(Pd_o))=abs((Pd_o_int-Pd_o))./Pd_o*100;
%     Error(green,1:length(Pd_g))=abs((Pd_g_int-Pd_g))./Pd_g*100;
%     Error(red,1:length(Pd_r))=abs((Pd_r_int-Pd_r))./Pd_r*100;
%     
    
    %%%loop
    for j=blue:red
        
        endline=find(SPDE_data(j,:)==0,1)
        if isempty(endline)
            endline=length(SPDE_data(j,:))
        end
        for i=1:endline;
%           
            final_list(i,2*j-1)=SPDE_data(j,i);
            final_list(i,2*j)=rms(stand_error(j,1:i));            
       end
   
    end
     plot(final_list(:,3),final_list(:,4),'m--',final_list(:,5),final_list(:,6),'g--',final_list(:,7),...
         final_list(:,8),'r--')
%     xlswrite('allDeviation.xlsx',final_list);
%     %%%%%%%%%%%%%%%
% semilogy(SPDE_simu(olive,:),Pd_simu(olive,:),'m',SPDE_simu(green,:),Pd_simu(green,:),'g',SPDE_simu(red,:),Pd_simu(red,:),'r');
% xlabel('Single-Photon Detection Efficiency %');
% ylabel('Dark Count Probability,Pd'); 
% legend('-30 degree','0 degree','20 degree','Location','southeast');
% hold on
% 
%     semilogy(SPDE_o,Pd_o,'m--',SPDE_g,Pd_g,'g--',SPDE_r,Pd_r,'r--')
% 
% figure;
% plot(SPDE_o,abs((Pd_o-Pd_o_int)./Pd_o*100),'m--',SPDE_b,abs((Pd_b-Pd_b_int)./Pd_b*100),'b--',...
%     SPDE_g,abs((Pd_g-Pd_g_int)./Pd_g*100),'g--',SPDE_r,abs((Pd_r-Pd_r_int)./Pd_r*100),'r--');
