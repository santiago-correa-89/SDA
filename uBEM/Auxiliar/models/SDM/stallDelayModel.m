function Cl = stallDelayModel(Clinv, Cl2d, thetaTwist, c, r)
% P. K. Chaviaropoulos and M. O. L. Hansen (2000) Stall Delay Model

thetaTwist = deg2rad(thetaTwist);
fCL        = 2.2*((c/r)^(1.3))*cos(thetaTwist)^4;
        
% Stall delay correction
Cl = Cl2d + fCL*( Clinv - Cl2d );

end