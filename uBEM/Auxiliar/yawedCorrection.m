function [Wx, Wz] = yawedCorrection(wx, wz, chi, theta_w, theta_yaw, r, R)

if theta_yaw < 0
    theta_w0 = pi/2;
elseif theta_yaw >= 0
    theta_w0 = -pi/2;
end

dum    = 1 + (r/R)*tan(chi/2)*cos(theta_w - theta_w0) ; 

Wx = wx*dum;
Wz = wz*dum;

end