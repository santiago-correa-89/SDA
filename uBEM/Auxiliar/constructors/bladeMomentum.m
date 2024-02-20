function M = bladeMomentum(r, Fld, R34)

R = length(r);
Raux  = zeros(3,R); 
M = 0; % Inicializo Momento de la pala

for j = 1:R
    Fld(:,j) = R34*Fld(:,j);
    Raux(:,j) = R34*[r(j); 0; 0];
end

for i = 1:(length(r) - 1)
    clear A; clear B;
    A   =  (Fld(2, i+1) - Fld(2, i)) / (Raux(1, i+1) - Raux(1, i));
    B   =  (Fld(2, i)*Raux(1, i+1) - Fld(2, i+1)*Raux(1, i)) / (Raux(1, i+1) - Raux(1, i));
    M   =  M + ( (1/3)*A*( Raux(1, i+1)^3 - Raux(1, i)^3 ) + (1/2)*B*( Raux(1, i+1)^2 - Raux(1, i)^2) );
end

% for i = 1:(length(r) - 1)
%     A   =  (Fld(2, i+1) - Fld(2, i)) / (r(i+1) - r(i));
%     B   =  (Fld(2, i)*r(i+1) - Fld(2, i+1)*r(i)) / (r(i+1) - r(i));
%     M   =  M + ( (1/3)*A*( r(i+1)^3 - r(i)^3 ) + (1/2)*B*( r(i+1)^2 - r(i)^2) );
% end

end