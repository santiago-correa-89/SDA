function [ Rcone ] = coneMatrix(thetaCone)

Rcone = [  cos(thetaCone)     sin(thetaCone)  0 ;...
           -sin(thetaCone)    cos(thetaCone)  0 ;...
                 0                   0        1 ];
           
end