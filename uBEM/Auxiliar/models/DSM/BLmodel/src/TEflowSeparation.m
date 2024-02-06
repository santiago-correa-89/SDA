function [ CnFn, CcFn, fn2primen, DFn ] = TEflowSeparation( sigma1, CcPn, CnCn, CnNCn, deltaS, ...
                                                      fnprimen, fnprimen1, DFn1, Tf0, eta) 

%% Viscouse component
Tf        = Tf0/sigma1 ;
DFn       = DFn1*exp(-(deltaS)/Tf) + (fnprimen - fnprimen1)*exp(-(deltaS)/(2*Tf));
fn2primen = fnprimen - DFn ;

%% Trailing edge flow separation components
CnFn  = CnCn*( ( (1 + sqrt(fn2primen) ) / 2 )^2 ) + CnNCn ;

CcFn  = CcPn*eta*( sqrt(fn2primen) - 0.2 );

end