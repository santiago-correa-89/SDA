function [ fnprimen, AoAfn, DPn, CnPrimen ] = effectiveFunction( Cn, CnPn, CnPn1, AoA0, AoA1, ...
                                                            AoA2, DPn1, deltaS, Tp, CnSlope, beta, ...
                                                            boolLeishmann )

S1   = 1; S2   = 1; S3    = 1; S4   = 1;

%% Effective angle of attack
DPn      = DPn1*exp( -deltaS/Tp )  +  (CnPn - CnPn1)*exp(-deltaS/(2*Tp));
CnPrimen = CnPn - DPn;
AoAfn    = ( CnPrimen / (CnSlope/beta) ) - AoA0 ;

%% Effective separation point
if ~boolLeishmann
    if AoA0 < AoAfn && AoAfn < AoA1
        fnprimen = 1 - 0.3*exp( (AoAfn - AoA1)/S1 ) ;
    elseif Aoa2 < AoAfn && AoAfn < AoA0
        fnprimen = 1 - 0.3*exp( (AoA2 - AoAfn)/S3 ) ;
    elseif AoAfn > AoA1
        fnprimen = 0.04 + 0.66*exp( (AoAf1 - AoAfn)/S2 ) ;
    elseif AoAfn < AoA2 
        fnprimen = 0.04 + 0.66*exp( (AoAfn - AoA2)/S4 ) ;
    end
else
    fnprimen  = max( min(( 2*sqrt( Cn/( (CnSlope/beta)*( AoAfn + AoA0) ) ) - 1 )^2 ,1 ), 0) ;
end