function [ bladePitchCnt, lastTimePC, speedError, integError ] = initPitchControl( nGen, omega, bladePitch, integError, lastTimePC, t1, V0 )      

bladePitchRad = deg2rad(bladePitch) ;

%% Control Parameter
%thetaK        = 6.302336    ; % Theta K  
thetaK        = 0.1099965         ; % Theta K radians
omegaN        = 0.6               ; % Natural frequency parameter of control
omegaSP       = 12.1*0.104719755  ; % angular velocity at TSR
Igen          = 534.116           ; % Generator inertia
Irot          = 115926            ; % Rotor inertia
Iblade        = 11776047          ; % Blade inertia in the root 
I             = nGen*nGen*Igen + Irot + 3*Iblade ; % Total inertia

dPdtheta = -25.52e6           ;
KI = I * omegaSP * omegaN * omegaN /nGen / (-dPdtheta) ;
GK = 1/ (1 + (bladePitchRad(1,1) / thetaK ) )                  ;
for i =1:3
    bladePitchCnt(1,i) = bladePitchRad(i)                  ;
    speedError    = nGen*omega - omegaSP           ;
    integError    = integError + bladePitchCnt(1,i)/( GK*KI )        ;
    bladePitchCnt(1,i) = rad2deg(bladePitchRad(1,i))         ;
end

end