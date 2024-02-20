function [ forces ] =  uBEMnormtanForces( bemForces, theta )

R = [1       0             0      ; ...  
     0   cos(theta)   sin(theta)  ; ...
     0   sin(theta)   -cos(theta)  ];

forces = R*bemForces;