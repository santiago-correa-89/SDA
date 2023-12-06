function T = thrustFactor(r, Pn)

T = 0; % Initializa blade thrust

for i = 1:(length(r) - 1)
    clear A; clear B; clear Ti;
    A   =  (Pn(i+1) - Pn(i)) / (r(i+1) - r(i));
    B   =  (Pn(i)*r(i+1) - Pn(i+1)*r(i)) / (r(i+1) - r(i));
    Ti  = ( (1/2)*A*( r(i+1)^2 - r(i)^2 ) + B*( r(i+1) - r(i) ) );
    T   =  T + Ti;
end
end