function [ genOmega, lastTimeGC, thetaPitch, integError, speedError, lastTimePC ] = initControlVariables(nblade, tIter)

% Torque control values
lastTimeGC       = zeros(tIter, 1)     ;
genOmega         = zeros(tIter, 1)     ;
% Pitch control values
integError       = zeros(tIter, 1)     ;
speedError       = zeros(tIter, 1)     ;
lastTimePC       = zeros(tIter, 1)     ;
thetaPitch       = zeros(tIter, nblade) ;

end