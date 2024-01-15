function [ CnPn, CcPn, CnCn, CnNCn, AoAEffn, Xn, Yn, Dn ]= attachedFlow( AoAn, deltaAoAn, KAoAn, KAoAn1, ...
                                                                    AoA0, CnSlope, deltaT, Dn1, M, beta, deltaS, ...
                                                                    Xn1, Yn1, Ti, A1, A2, b1, b2)

% Circulatory component
Xn        = Xn1*exp(-b1*( beta^2 )*deltaS) + A1*deltaAoAn*exp(-b1*( beta^2 )*deltaS/2); % Indical X function
Yn        = Yn1*exp(-b2*( beta^2 )*deltaS) + A2*deltaAoAn*exp(-b2*( beta^2 )*deltaS/2); % Indical Y function
AoAEffn   = AoAn - Xn - Yn;                                                    % Effective angle of attack
CnCn      = ( CnSlope/beta )*(AoAEffn + AoA0);                                                   % Ciculatori AoA component

% Non circulatory component
Ka        = ( ( 1 - M ) + CnSlope*beta*(M^2)*(A1*b1 + A2*b2) )^(-1);                 % Non circulatory AoA constante
TAoA      = 0.75*Ka*Ti ;                                                          % Non circulatory time constante
Dn        = Dn1*exp( -deltaT/TAoA ) + ( KAoAn - KAoAn1 )*exp( -deltaT/(2*TAoA) ); % Non circulatory defciency function
CnNCn     = ( (4*TAoA)/M )*( KAoAn - Dn );                                        % Non circulatory normal force

% Total normal force under attached conditions
CnPn   = CnCn + CnNCn;

% Total chord force under attached conditions
CcPn   = CnCn*tan( AoAEffn + AoA0);