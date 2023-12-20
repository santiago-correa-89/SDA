function [ Cl, Cd, Cm ] = airfoilCoef( alpha, datos)
%% Función que permite determinar valores de sustentación y arrastre
%% Search for airfoil data

%%Load file data
alphaData = datos(:,1);
ClData    = datos(:,2);
CdData    = datos(:,3);
CmData    = datos(:,4);

Cl        = interp1(alphaData, ClData, rad2deg(alpha), 'linear', 'extrap');
Cd        = interp1(alphaData, CdData, rad2deg(alpha), 'linear', 'extrap');
Cm        = interp1(alphaData, CmData, rad2deg(alpha), 'linear', 'extrap');

end