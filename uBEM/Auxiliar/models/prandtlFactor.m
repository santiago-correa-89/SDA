function F = prandtlFactor(r, R, Rhub, alpha)
      
    f_tip = 3*(R - r)/(2*r*abs(sin(alpha)));
    F_tip = (2/pi)*acos(exp(-f_tip));
            
    f_hub = 1; %3*( r - Rhub )/(2*r*abs(sin(alpha)));
    F_hub = 1; %2/pi)*acos(exp(-f_hub));

    F     = F_tip * F_hub;
end