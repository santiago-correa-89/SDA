function [ bladePitchOut, lastTimePC, speedError, integError ] = pitchControl( nGen, rotOmega, bladePitch, integError, lastTimePC, t1p, V0 )      

bladePitchOut = zeros(1,3) ;
bladePitchRad = deg2rad(bladePitch) ;

%% Control Parameter
minPit        = 0           ; % Minimum pitch setting
maxPit        = 1.570796    ; % Maximum pitch setting
%thetaK        = 6.302336    ; % Theta K  
thetaK        = 0.1099965        ;
seda          = 0.7              ; % Damping control paramater
omegaN        = 0.6                   ; % Natural frequency parameter of control
maxPitchRate  = 8*pi/180              ; % Theta pitch change rate
omegaSP       = 12.1*0.104719755      ; % angular velocity at TSR
Igen          = 534.116               ; % Generator inertia
Irot          = 115926                ; % Rotor inertia
Iblade        = 11776047               ; % Blade inertia in the root 
I             = nGen*nGen*Igen + Irot + 3*Iblade ; % Total inertia

dPdtheta   = -25.52e6;
GK = 1/ (1 + (bladePitchRad(1,1) / thetaK ) ) ;
KP = ( 2 * I * omegaSP * seda * omegaN) / nGen / (-dPdtheta)   ;
KI = I * omegaSP * omegaN * omegaN /nGen / (-dPdtheta) ;

deltaTime     = t1p - lastTimePC ; % Delta between control check

    genOmega   = nGen*rotOmega ;
    %cornerFreq = 1.570796                                  ; % Frequency parameter of filter
    %alpha      = exp( ( deltaTime )*cornerFreq )           ; % Constant parameter of filter
    genOmegaF  = genOmega                                  ; % filter gen speed

    %% Set control parameters
    speedError = - nGen*omegaSP + genOmegaF ;              % speed error to compute PID control
    integError = integError + speedError*deltaTime; % acumulated speed error to compute PID control

    %% Saturate the integral error term
    integError = min( max( integError, minPit / (GK * KI) ) , maxPit / (GK * KI) );

    %% Compute control output
    pitchConP = GK*KP*speedError ;
    pitchConI = GK*KI*integError ;

    bladePitchCnt  = pitchConP + pitchConI                         ;
    bladePitchCnt  = min( max( bladePitchCnt, minPit), maxPit)     ; % Saturate the overall command using the pitch angle
    
    for i = 1:3
        pitchRate = ( bladePitchCnt - bladePitchRad(1,i) )/deltaTime            ;
        pitchRate = min ( max(pitchRate, -maxPitchRate), maxPitchRate)        ;
    
        bladePitchOut(1,i) = bladePitchRad(1,i) + pitchRate*deltaTime          ;  % compute blade pitch variation
        bladePitchOut(1,i) = rad2deg( bladePitchOut(1,i))                     ;
    end
    lastTimePC = t1p ;

end


