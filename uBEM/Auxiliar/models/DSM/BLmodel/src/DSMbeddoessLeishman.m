function [BLcoefs] = DSMbeddoessLeishman(  AoAn, Vrel, deltaT, c,  BLcoefs, ...
                                          AoA0,  Cd0, CnSlope, Clst, Cdst, Cn1, Cn2, AoA1, AoA2, T, boolLeishmann )

% AoAn     - current angle of attack at step n
% Vrel     - current relative velocity of the section
% deltaT   - Time ste
% c        - Chord of the section
% BLcoefs  - Beddoess and Leishman coef of previous time step (n - 1)
% AoA0     - Zero lift angle of attack of the section

AoAn  = deg2rad( AoAn )  ; AoA0  = deg2rad( AoA0 )  ; AoA1  = deg2rad( AoA1 )  ; AoA2  = deg2rad( AoA2 )  ;

% Initialize parameters
a        = 340;                  % Sound Velocity
deltaS   = (Vrel*deltaT)/(c/2);  % Función de distancia adimensional
M        = Vrel/a;               % Numero de Mach
beta     = sqrt(1 - M^2);        % Constante beta del número de Mach
CnCSlope = CnSlope/beta       ;  % Slope of the circulatory normal componente

% Parámetros y constantes de los de los modelo
A1 = 0.3 ;  A2  = 0.7 ;   b1  = 0.14;  b2  = 0.53;
Tp = T(1);  Tf0 = T(2);   Tv0 = T(3);  Tvl = T(4);      Ti = c/a;   
St = 0.19;  eta = 0.95;
sedaLP          = 0.5;                         % Low pass filter frequency cutoff -3DB      
Clp             = exp(-2*pi*deltaT*sedaLP);    % Low pass filter constant

%%Initial conditions
% Data of flux at previous step
Cl1         = BLcoefs(1)         ; % Cl coef estimation with LBM at previous time step 
Cd1         = BLcoefs(2)         ; % Cd coef estimation with LBM at previous time step
CNn1        = BLcoefs(3)         ; % Cn coef estimation with LBM at previous time step
CCn1        = BLcoefs(4)         ; % Cc coef estimation with LBM at previous time step
AoAEffn1    = BLcoefs(5)         ; % Effective angle of attack of previous time step
AoAn1       = BLcoefs(6)         ; % Filtered Angle of attack of previous time step
qn1         = BLcoefs(7)         ; % Filtered non dimensional pitching rate of previous time step
deltaAoAn1  = BLcoefs(8)         ; % Filtered Angle of attack variation at previous time step
KAoAn1      = BLcoefs(9)         ; % Filtered Discretd derivative of angle of attack

% Attached flow data at previous step
Xn1         = BLcoefs(10)        ; % Xn circulatory component at previous time step
Yn1         = BLcoefs(11)        ; % Yn circulatory component at previous time step
Dn1         = BLcoefs(12)        ; % Dn deficeny function at previous time step

% Trailing edge flow separation parameters of previous step
DPn1        = BLcoefs(13)        ; % Dpn deficeny function of attached normal force at previous time step
DFn1        = BLcoefs(14)        ; % Dfn of deficeny function of separation point at previous time step
CnPn1       = BLcoefs(15)        ; % Normal force due to attached flow at previous time step 
fnprimen1   = BLcoefs(16)        ; % Effective separation point at previous time step

% Flow separated data at leading edge due to vortical separation
fn2primen1  = BLcoefs(17)        ; % Effective separation point
AoAfn1      = BLcoefs(18)        ; % Effective incidence angle at previous time step
CnVn1       = BLcoefs(19)        ; % Additional lift force due to vortex shedding at previous time step
CVn1        = BLcoefs(20)        ; % Accumulated vorticity at LE contribution
tauVn1      = BLcoefs(21)        ; % Non dimensional time variable 
CnFn1       = BLcoefs(22)        ;


% Low pass filter AoA, q Kq
%AoAn       = Clp*AoAn1 + (1 - Clp)*AoAn     ; % Filtered angle of attack
qn         = (AoAn - AoAn1)*c/(Vrel*deltaT) ; % Non dimensional pitching rate
%qn         = Clp*qn1 + (1 - Clp)*qn         ; % Fltered Non dimensional pitching rate
deltaAoAn  = ( AoAn - AoAn1 )               ; % Filtered angle of attack variation
KAoAn      = qn*Vrel/c                      ; % Filtered angle of attack derivative
deltaAoA0  = AoAn - AoA0                    ; % Filtered angle of attack variation

%% Static normal and chord coefs
Cn  = Clst*cos(AoAn)  +  (Cdst - Cd0)*sin(AoAn) ;
Cc  = Clst*sin(AoAn)  -  (Cdst - Cd0)*cos(AoAn) ;

%% Attached Flow Module
[ CnPn, CcPn, CnCn, CnNCn, AoAEffn, Xn, Yn, Dn ] = attachedFlow( AoAn, deltaAoAn, KAoAn, KAoAn1, ...
                                                                    AoA0, CnSlope, deltaT, Dn1, M, beta, deltaS, ...
                                                                    Xn1, Yn1, Ti, A1, A2, b1, b2);

[ fnprimen, AoAfn, DPn, CnPrimen ]  = effectiveFunction( Cn, CnPn, CnPn1, AoA0, AoA1, ...
                                                AoA2, DPn1, deltaS, Tp, CnSlope, beta, ...
                                                boolLeishmann);

%% Trailing edge flow separation
sigma1 = 1; % Default value
[ CnFn, CcFn, fn2primen, DFn ] = TEflowSeparation( sigma1, CcPn, CnCn, CnNCn, deltaS, ...
                                                      fnprimen, fnprimen1, DFn1, Tf0, eta);

%% Dynamic stall vortex advection
% if fn2primen > 1
%     fn2primen = 1 ;
% elseif fn2primen < 0
%     fn2primen = 0 ;
% end
Tsh   = (2*( 1 - fn2primen ))/St;        % Shedding frequency
Tvl   = min( max(Tvl, 6), 13);           % Min Max range of Fast8 documentation

KNn   = (1 + sqrt( fn2primen ) )/ 2;     % 
CVn   = CnCn*( (1 - KNn)^2 );            % Contribution of normal force due to voricity accumulation

% Vortex non dimensional displement 
tauVn   = tauVn1 + ( 2*Vrel*deltaT/c );

%% Flow conditions flags
% LE condition
if CnPrimen > Cn1 && AoAn >= AoA0 || CnPrimen < Cn2 && AoAn < AoA0
    LESF = true; % LE separation take place
else
    LESF = false; % Reattachment can occur
end
% Vortex condition
if 0 < tauVn && tauVn <= 2*Tvl
    VRTX = true;  % Vortex advection in progress
else
    VRTX = false; % Vortex in the wake
end

if tauVn > Tvl + Tsh && LESF
    tauVn = 0;
end

%% Tf modifications
if fn2primen < fn2primen1    % TE separation take place
    TESF = true; 
    if KAoAn*deltaAoA0 < 0   % accelerate separation point movement
        sigma1 = 2;
    elseif ~LESF             % default value, LE separation could take place
        sigma1 = 1 ;
    elseif LESF              % Leading edge flow separation has begun resulting in
        sigma1 = 1.75;       % accelerated trailing edge separation point movement.
    elseif fn2primen1 <= 0.7 % accelerate separation point movement if separation is occurring) 
        sigma1 = 2;
    end
else 
    TESF = false ;          % Flow reattachment
    if ~LESF 
        sigma1 = 0.5;       % Default: slow down reattachment
    elseif 0 < tauVn && tauVn <= Tvl
        sigma1 = 0.25;      % No flow reattachment if vortex shedding is in progress
    elseif KAoAn*deltaAoA0 > 0
        sigma1 = 0.75;
    end
end

%% Update Trailing edge flow separation
[ CnFn, CcFn, fn2primen, DFn ] = TEflowSeparation( sigma1, CcPn, CnCn, CnNCn, deltaS, ...
                                                      fnprimen, fnprimen1, DFn1, Tf0, eta);

%% TV modifications
sigma3 = 1;                            % Initalize default value sigma 3
if Tvl < tauVn && tauVn <= 2*Tvl       % Postshedding
    sigma3 = 3; 
elseif ~TESF                           % Accelerate vortex lift decay
    sigma3 = 4;
elseif VRTX && 0 <= tauVn && tauVn <= Tvl  % Vortex shedding (Default)
    sigma3 = 1;
    if KAoAn*deltaAoA0 < 0     % Accelerate vortex lift decay
        sigma3 = 2;
    end
elseif KAoAn*deltaAoA0 < 0             % vortex lift must decay fast
    sigma3 = 4;
    if ~TESF   % flow is reattaching and angle of attack is decreasing (Default) 
        sigma3 = 1;
    end
end

%% Vortex advection model
[ CnVn ] = LEvortexAdvection(sigma3, tauVn, CnVn1, CVn, CVn1, ...
                                        deltaS, Tvl, Tv0, LESF);

%% Determination of Cn, Cc, Cl and Cd
CNn  = CnFn  +  CnVn ;

CCn = eta*(CnSlope/beta)*(AoAEffn)^2*sqrt(fn2primen);

%CCn = CcFn + CnVn*tan(AoAEffn)*( 1 - tauVn/Tvl ) ;

% if CnPrimen <= Cn2
%     CCn   = CcPn*eta*( sqrt(fn2primen) )*cos( AoAEffn ) ;
% else
%     phi   = 2*(CnPrimen - Cn2) + (fn2primen - fnprimen) ;
%     CCn   = CcPn*eta*( sqrt(fn2primen) )*( fn2primen^(phi) )*cos( AoAEffn ) ;
% end

% Total (Ejes viento)
Cln = CNn*cos(AoAn)  +  CCn*sin(AoAn);
Cdn = CNn*sin(AoAn)  -  CCn*cos(AoAn) + Cd0;

%% Update BL coefs 

BLcoefs(1)  = Cln       ;
BLcoefs(2)  = Cdn       ;
BLcoefs(3)  = CNn       ;
BLcoefs(4)  = CCn       ;
BLcoefs(5)  = AoAEffn   ;  % Effective angle of attack of previous time step
BLcoefs(6)  = AoAn      ;  % Filtered Angle of attack of previous time step
BLcoefs(7)  = qn        ;  % Filtered non dimensional pitching rate of previous time step
BLcoefs(8)  = deltaAoAn ;  % Angle of attack variation at previous time step
BLcoefs(9)  = KAoAn     ;  % Discretd derivative of angle of attack
BLcoefs(10) = Xn        ;  % Constante Componente Circulatoria
BLcoefs(11) = Yn        ;  % Constante Componente Circulatoria
BLcoefs(12) = Dn        ;  % Constante Componente no Circulatoria
BLcoefs(13) = DPn       ;
BLcoefs(14) = DFn       ;
BLcoefs(15) = CnPn      ;   
BLcoefs(16) = fnprimen  ;
BLcoefs(17) = fn2primen ;
BLcoefs(18) = AoAfn     ;
BLcoefs(19) = CnVn      ;
BLcoefs(20) = CVn       ;
BLcoefs(21) = tauVn     ;
BLcoefs(22) = CnFn      ;

end