function [Ct] = thrustFactor(r, Pn, rho, V0, Rrot, N)

T = 0; % Inicializo Momento de la pala

for i = 1:(length(r) - 1)
    clear A; clear B; clear Ti;
    A   =  (Pn(i+1) - Pn(i)) / (r(i+1) - r(i));
    B   =  (Pn(i)*r(i+1) - Pn(i+1)*r(i)) / (r(i+1) - r(i));
    Ti  = ( (1/2)*A*( r(i+1)^2 - r(i)^2 ) + B*( r(i+1) - r(i) ) );
    T   =  T + Ti;
end

Ct = N*T/((1/2)*rho*(V0^2)*pi*(Rrot^2));

end