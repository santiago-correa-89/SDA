function [M, Cp] = powerFactor(r, Pt, omega, rho, V0, Rrot, N)

M = 0; % Inicializo Momento de la pala

for i = 1:(length(r) - 1)
    A   =  (Pt(i+1) - Pt(i)) / (r(i+1) - r(i));
    B   =  (Pt(i)*r(i+1) - Pt(i+1)*r(i)) / (r(i+1) - r(i));
    M   =  M + ( (1/3)*A*( r(i+1)^3 - r(i)^3 ) + (1/2)*B*( r(i+1)^2 - r(i)^2) );
end

M  = N*M;
Cp = omega*M/((1/2)*rho*(V0^3)*pi*(Rrot^2));

end