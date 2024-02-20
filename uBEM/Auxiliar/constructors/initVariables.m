function [ Vblade, rG, V0, thetaWing, omega, wqs, wint, wake, force, clift, cdrag, AoA, flowAngle, fs ] = initVariables(ni, nblade, tIter)

Vblade          = zeros( 3, tIter )             ;
rG              = zeros( 3, ni, nblade,  tIter) ;
V0              = zeros( 3, ni, nblade,  tIter) ;
wqs             = zeros( 3, ni, nblade,  tIter) ;
wint            = zeros( 3, ni, nblade,  tIter) ;
wake            = zeros( 3, ni, nblade,  tIter) ;
thetaWing       = zeros( tIter, 3)              ;
omega           = zeros( tIter, 1)              ;
force           = zeros( 3, ni, nblade,  tIter) ;
clift           = zeros( ni, nblade,  tIter)    ;
cdrag           = zeros( ni, nblade,  tIter)    ;
AoA             = zeros( ni, nblade,  tIter)    ;
flowAngle       = zeros( ni, nblade,  tIter)    ;
fs              = zeros( ni, nblade,  tIter)    ; % if oye model is turned on