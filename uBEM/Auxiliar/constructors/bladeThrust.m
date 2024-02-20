function T = bladeThrust(r, Fld, Rcone)

T = 0; % Initializa blade thrust

for i = 1:(length(r) - 1)
    clear A; clear B; clear Ti;
    Fth1 =  Rcone*Fld(:,i)   ;
    Fth2 =  Rcone*Fld(:,i+1) ;
    A    =  (Fth2(3) - Fth1(3)) / (r(i+1) - r(i));
    B    =  (Fth1(3)*r(i+1) - Fth2(3)*r(i)) / (r(i+1) - r(i));
    Ti   = ( (1/2)*A*( r(i+1)^2 - r(i)^2 ) + B*( r(i+1) - r(i) ) );
    T    =  T + Ti;
end
end