function [ CnVn ] = LEvortexAdvection(sigma3, tauVn, CnVn1, CVn, CVn1, ...
                                        deltaS, Tvl, Tv0, LESF)

Tv   = Tv0 / sigma3 ;

if 0 <= tauVn && tauVn <= Tvl
    CnVn   = CnVn1*exp(-deltaS/Tv) + (CVn - CVn1)*exp(-deltaS/(2*Tv));

elseif tauVn > Tvl || ~LESF
    sigma3 = 2;
    Tv     = Tv0/sigma3 ;
    CnVn   = CnVn1*exp( -deltaS/Tv ) ;
end
