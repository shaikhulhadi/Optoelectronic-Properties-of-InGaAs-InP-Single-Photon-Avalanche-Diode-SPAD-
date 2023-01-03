clear all;
close all;
clc;
%Temperature=30;
Vapplied=10:.2:67.1;
Applied_Length=length(Vapplied);
Plank = 6.63e-34;
Light=3e8;
Electron=1.6e-19;
Lambda = 1550e-9;
Input_Power = 196e-9;
Vbreakdown=67.1;
Vpunchthrough=35;%%Chosen arbitrarily
n=3.5; 
R = .65*(Electron*Lambda)/(Plank*Light);
Ipho = R*Input_Power;

red_20c=csvread('red.csv');

Iph_expe=red_20c(:,2)';
Voltage_expe=red_20c(:,1)';
% M= 1./(1-(Vapplied./Vbreakdown).^n

for i=1:length(Vapplied)
    
    T(i)=Vapplied(i)/Vbreakdown;
    
    
    
    if (Vapplied(i)<=62)
        
        M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)));

    end      
             
     if (Vapplied(i)>62) && (Vapplied(i)<=66)

             M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)));

        
     end
     
      if (Vapplied(i)>66)
          
              M(i)=1+((T(i))^n)+(((T(i))^(2*n)))+(((T(i))^(3*n)))+(((T(i))^(4*n)))+(((T(i))^(5*n)))+...
            (((T(i))^(6*n)))+(((T(i))^(7*n)))+(((T(i))^(8*n)))+(((T(i))^(9*n)))+(((T(i))^(10*n)))+...
            (((T(i))^(11*n)));
     
     
      end
     
    end

    

 
I_out = M*Ipho; 

plot(Vapplied,I_out,'k',Voltage_expe,Iph_expe,'r');