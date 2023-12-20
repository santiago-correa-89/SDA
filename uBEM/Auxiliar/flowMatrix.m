function [ Rangle ] = flowMatrix(thetaFlow)

Rangle = [  cos(thetaFlow)     0     sin(thetaFlow) ;...
                 0             1           0          ;...
            -sin(thetaFlow)    0      cos(thetaFlow)  ];
           
end