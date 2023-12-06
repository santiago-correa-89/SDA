function chi = skewAngle(V0, W, a34, a13)

Wyaw   = a34'*W ;
V0l3   = a13*V0 ; 

% Definir el vector n en el sistema de coordenadas del rotor YAwed
n_yaw  = [1; 0; 0];  % Por ejemplo, n_inclinado está a lo largo del eje x

Vyawed = V0l3 + n_yaw*dot(n_yaw,Wyaw) ;

chi    = acos(dot(Vyawed, n_yaw)/norm(Vyawed)) ;

