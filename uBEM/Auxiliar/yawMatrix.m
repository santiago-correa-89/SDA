function [ Ryaw ] = yawMatrix(thetaYaw)

Ryaw = [  cos(thetaYaw)     0    sin(thetaYaw) ;  ...
               0            1          0       ;  ...
          -sin(thetaYaw)    0    cos(thetaYaw)  ];
end