function [Cl, Cd] = airfoilCoef( alpha, datos)
%% Función que permite determinar valores de sustentación y arrastre
%% Search for airfoil data

%%Load file data
alphaData = datos.data(:,1);
ClData    = datos.data(:,2);
CdData    = datos.data(:,3);

Cl        = interp1(alphaData, ClData, rad2deg(alpha), 'linear', 'extrap');
Cd        = interp1(alphaData, CdData, rad2deg(alpha), 'linear', 'extrap');

end