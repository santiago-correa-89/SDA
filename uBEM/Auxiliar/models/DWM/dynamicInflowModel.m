function [Wake, WakeInt] = dynamicInflowModel(WakeOld, WakeQSOld, WakeIntOld, WakeQS, normVrel, r, R, a, R0, deltaT)

% Dynamic Wake Model (Snel and Schepers, 1995). Oye Model

% Express wake variables into current coordinate system
WakeOld     =  R0*WakeOld    ;
WakeQSOld   =  R0*WakeQSOld  ;
WakeIntOld  =  R0*WakeIntOld ;

k = 0.6; 

Tau1 = ( 1.1 / ( 1 - 1.3*a ) ) * R / normVrel;
Tau2 = ( 0.39 - ( 0.26*(r/R)^2 ) )*Tau1;

% Calculate right hand of first order differential equation using backward
% difference
H = WakeQS + k*Tau1*(WakeQS - WakeQSOld) / deltaT;
              
% Solve intermedial Induced Velocity solving first differential equation
WakeInt = H + ( WakeIntOld - H ) * exp(-deltaT/Tau1);
              
% Solve analitically second first order differential equation        
Wake = WakeInt + ( WakeOld - WakeInt ) * exp(-deltaT/Tau2);
          
end
