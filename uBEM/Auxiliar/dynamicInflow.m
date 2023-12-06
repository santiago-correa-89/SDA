function [Wx, Wz, Wxint, Wzint] = dynamicInflow(Wxold, Wxqs, Wxqsold, Wxintold, Wzold, Wzqs, Wzqsold, Wzintold, Vrel, r, R, a, deltaT)

% Dynamic Wake Model (Snel and Schepers, 1995). Oye Model

k = 0.6; 

Tau1 = ( 1.1 / ( 1 - 1.3*a ) ) * R / sqrt(Vrel);
Tau2 = ( 0.39 - ( 0.26*(r/R)^2 ) )*Tau1;

% Calculate right hand of first order differential equation using backward
% difference
Hx = Wxqs + k*Tau1*(Wxqs - Wxqsold) / deltaT;
Hz = Wzqs + k*Tau1*(Wzqs - Wzqsold) / deltaT;
              
% Solve intermedial Induced Velocity solving first differential equation
Wxint = Hx + ( Wxintold - Hx ) * exp(-deltaT/Tau1);
Wzint = Hz + ( Wzintold - Hz ) * exp(-deltaT/Tau1);
              
% Solve analitically second first order differential equation        
Wx = Wxint + ( Wxold - Wxint ) * exp(-deltaT/Tau2);
Wz = Wzint + ( Wzold - Wzint ) * exp(-deltaT/Tau2);
          
end
