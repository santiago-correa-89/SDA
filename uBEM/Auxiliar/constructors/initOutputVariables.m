function [ aeroPower, genPower, blThrust, thrust, genTrq, aeroTrq, rotTrq, CL ] = initOutputVariables(tIter, nblade)

% Outputs
aeroPower     = zeros(tIter, 1) ;
genPower      = zeros(tIter, 1) ;
genTrq        = zeros(tIter, 1) ;
rotTrq        = zeros(tIter, 1) ;
thrust        = zeros(tIter, 1) ;
CL            = zeros(tIter, 1) ;
aeroTrq       = zeros(tIter, nblade) ;
blThrust      = zeros(tIter, nblade) ; 

end