%% Método BEM NREL 5 MW
clear all; 
close all; 
clc;

%% Párametros del problema
% theta_t = ángulo de giro (grados)
% theta_p = ángulo de pitch (grados)
% B       = Número de palas
% c       = cuerda
% r       = coordenada radial
% R       = radio de la pala
% lambda  = tip speed ratio
% Cl      = lift coef
% Cd      = drag coef
% a       = axial indutcion coef
% a_prime = induced velocity coef
% phi     = angulo de velocidad relativa
% err     = error de iteración
% k       = número de iteraciones
%Dhub     = Diametro del hub

%% Importar data de la pala
% path          = 'data/Tjaereborg' ; 
% file          = '/Tjaereborg_2.0MW_061m.txt';

addpath( genpath( [ pwd '/BEM']  ) );
addpath( genpath( [ pwd '/uBEM/Data'] ) );

path       = '/NREL_5MW' ;
file       = '/Polars/NRELOffshrBsline5MWbladeDataV1.txt';
file1      = '/Polars/Cylinder1.txt' ;
file2      = '/Polars/Cylinder2.txt' ;
file3      = '/Polars/DU21_A17.txt'  ;
file4      = '/Polars/DU25_A17.txt'  ;
file5      = '/Polars/DU30_A17.txt'  ;
file6      = '/Polars/DU35_A17.txt'  ;
file7      = '/Polars/DU40_A17.txt'  ;
file8      = '/Polars/NACA64_A17.txt';

polars = { file1, file2, file3, file4, file5, file6, file7, file8 } ;

% Abre el archivo en modo de lectura ('r' para lectura)
bladeData   = readtable(fullfile(path, file));
radius      = bladeData.RNodes    ;
twist       = bladeData.AeroTwst  ; % Angle in DEGRESS
chord       = bladeData.Chord     ;
sectionID   = bladeData.NFoil     ;
ni          = max(size(radius))   ;

% check airfoil data of the section using the ID 
[aoa, clift, cdrag, cmomt ] = readairfoildata_v3( polars, path ) ;

% Comprueba si el archivo se abrió correctamente
if isempty(bladeData)
    error('No se pudo abrir el archivo.');
end

%% Datos del problema
Rhub    = 3/2   ; % Radio del Hub
Rrot    = 126/2 ; % Radio del rotor

theta_p = deg2rad(-5:1:5)' ; % Vector de Ángulo de Twist de la las secciones en radianes

V0     = 11.4;              % m/s
lambda = (3:1:12)';         % Tip speed ratio
omega  = (lambda*V0/Rrot) ; % Velocidad angular
B      = 3;                 % Número de palas
rho    = 1.225;             % Densidad del aire 

%% Párametros de Iteración
% Factor de corrección por factor de inducción axial "a" elevado
% factorCorreccion = 0 - No aplica factor de corrección
% factorCorreccion = 1 - Factor de corrección de Glauert
% factorCorreccion = 2 - Factor de corrección de Wilson and Walker 1984 (Spera 1994)
factorCorreccion = 1;
% Tolerancia de iteración BEM
Tol = 1e-3;

%% Reservo memoria para parámetros
a        = zeros(length(radius), length(lambda), length(theta_p));
aprime   = zeros(length(radius), length(lambda), length(theta_p));
phi      = zeros(length(radius), length(lambda), length(theta_p));
alpha    = zeros(length(radius), length(lambda), length(theta_p));
Vrel     = zeros(length(radius), length(lambda), length(theta_p));
lambda_r = zeros(length(radius), length(lambda), length(theta_p));
Cl       = zeros(length(radius), length(lambda), length(theta_p));
Cd       = zeros(length(radius), length(lambda), length(theta_p));
Pn       = zeros(length(radius), length(lambda), length(theta_p));
Pt       = zeros(length(radius), length(lambda), length(theta_p));
kiter    = zeros(length(radius), length(lambda), length(theta_p));
err      = zeros(length(radius), length(lambda), length(theta_p));
M        = zeros(length(lambda), length(theta_p));
Cp       = zeros(length(lambda), length(theta_p));
Ct       = zeros(length(lambda), length(theta_p));

%% Resolución de método BEM
for g = 1:length(theta_p)
    for h = 1:length(lambda)
        clear Ke;
        for i = 1:length(radius)
    
            r_i               = radius(i);        % Radio de pala
            c_i               = chord(i);         % Cuerda de la sección
            theta_t_i         = twist(i);         % Angulo de Pitch
            t_i               = sectionID(i);     % Dependiendo del codigo el espesor de la función permite identificar el perfil, si tenemos el ID cargamos en t el Id y lo usamos para levantar las polares
            
            lambda_r(i, h, g) = lambda(h)*r_i/Rrot;    % Cociente Velocidad del aire y velocidad angular en la sección

            [ a(i, h, g) , aprime(i, h, g) , phi(i, h, g) , alpha(i, h, g), Cl(i, h, g), Cd(i, h, g), err(i, h, g) , kiter(i, h, g) ] = BEM_method(B, c_i, r_i, t_i, ...
                                                                                                        lambda(h), theta_t_i, theta_p(g), Rrot, Rhub, Tol, 1, 1, factorCorreccion, aoa, clift, cdrag, cmomt );
   
        end

        Vrel(:, h, g)  = sqrt( (V0*(1-a(:, h, g))).^2 + ((lambda_r(:,h , g)*V0.*(1+aprime(:, h, g))).^2 ) );
        
        Ke             = (1/2)*rho*(Vrel(:, h, g).^2).*chord(:);
        Pn(:, h, g)    = Ke.*(Cl(:, h, g).*cos(phi(:, h, g)) + Cd(:, h, g).*sin(phi(:, h, g)));
        Pt(:, h, g)    = Ke.*(Cl(:, h, g).*sin(phi(:, h, g)) - Cd(:, h, g).*cos(phi(:, h, g)));
    end
end

%% Variables of BEM to estimate performance of HAWT

for i = 1:length(theta_p)
    for j= 1:length(lambda)
        [ M(j,i) , Cp(j,i) ] = powerFactor(radius, Pt(:, j, i), omega(j), rho, V0, Rrot, B) ;
        [ Ct(j,i) ]          = thrustFactor(radius, Pn(:, j, i), rho, V0, Rrot, B)          ;
    end
end

% Crear una cuadrícula más densa para la interpolación
lambda_interp  = linspace(min(lambda), max(lambda), length(lambda)*10); % Ajusta el número de puntos según tus necesidades
theta_p_interp = linspace(min(rad2deg(theta_p)), max(rad2deg(theta_p)), length(theta_p)*10); % Ajusta el número de puntos según tus necesidades
[lambda_interp_grid, theta_p_interp_grid] = meshgrid(lambda_interp, theta_p_interp);

% Interpolar los valores de Cp en la cuadrícula densa
Cp_interp = griddata(rad2deg(theta_p), lambda , Cp, theta_p_interp_grid, lambda_interp_grid, 'cubic'); % Puedes ajustar el método de interpolación
Ct_interp = griddata(rad2deg(theta_p), lambda , Ct, theta_p_interp_grid, lambda_interp_grid, 'cubic'); % Puedes ajustar el método de interpolación

% Plot results
folderPathFigs = './figs/NREL_5MW/' ;
mkdir(folderPathFigs) ;

lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ; 

% Crear un gráfico de superficie suavizado con contourf
fig1 = figure(1);
labelTitle = 'Grafico de superficie Cp' ;
hold on; grid on
contourf(lambda_interp_grid, theta_p_interp_grid, Cp_interp, 20, '--', "ShowText",true,"LabelFormat", "%0.2f"); % Ajusta el número de niveles según tus necesidades
labx = xlabel('\lambda'); laby = ylabel('\theta_p (grados)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig1 = strcat(folderPathFigs, 'superficeCp.png') ;
print(fig1, namefig1,'-dpng') ;

% Crear un gráfico de superficie suavizado con contourf
fig2 = figure(2);
labelTitle = 'Grafico de superficie C_t' ;
hold on; grid on
contourf(lambda_interp_grid, theta_p_interp_grid, Ct_interp, 20, '--', "ShowText",true,"LabelFormat", "%0.2f"); % Ajusta el número de niveles según tus necesidades
labx = xlabel('\lambda'); laby = ylabel('\theta_p (grados)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig2 = strcat(folderPathFigs, 'superficeCt.png') ;
print(fig2, namefig2,'-dpng') ;


theta_p_labels = { '\theta_p = -5', '\theta_p = -4', '\theta_p = -3', '\theta_p = -2', '\theta_p = -1', '\theta_p = 0', '\theta_p = 1', '\theta_p = 2', '\theta_p = 3', '\theta_p = 4', '\theta_p = 5' };
lambda_labels  = {'\lambda = 7', '\lambda = 9', '\lambda = 11'};

fig3 = figure(3);
labelTitle = 'Grafico Cp-\lambda' ;
hold on; grid on
for l =1:length(theta_p)
    figure(3)
    aux = length(Cp_interp)/length(theta_p);
    plot(lambda_interp, Cp_interp((l-1)*aux+1,:))
    hold on
end
labx = xlabel('\lambda'); laby = ylabel('Cp');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
legend(theta_p_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
title(labelTitle)
namefig3 = strcat(folderPathFigs, 'CpLambda.png') ;
print(fig3, namefig3,'-dpng') ;

fig4 = figure(4);
labelTitle = 'Grafico Ct-\lambda' ;
hold on; grid on
for l =1:length(theta_p)
    figure(4)
    aux = length(Ct_interp)/length(theta_p);
    plot(lambda_interp, Ct_interp((l-1)*aux+1,:))
    hold on
end
labx = xlabel('\lambda'); laby = ylabel('Ct');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
legend(theta_p_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
title(labelTitle)
namefig4 = strcat(folderPathFigs, 'CtLambda.png') ;
print(fig4, namefig4,'-dpng') ;


for i = 1:length(theta_p)
    for h = 1:length(lambda)
        r_R = radius/Rrot;
        if theta_p(i) == 0
            if lambda(h) == 7 || lambda(h) == 9 || lambda(h) == 11
                
                fig5 = figure(5);
                subplot(3,1,1)
                hold on
                plot(r_R, a(:, h, i))
                xlabel('r/R')
                ylabel('a')
                subplot(3,1,2)
                hold on
                plot(r_R, aprime(:, h, i))
                xlabel('r/R')
                ylabel('a_prime')
                subplot(3,1,3)
                hold on
                plot(r_R, alpha(:, h, i))
                xlabel('r/R')
                ylabel('AoA (Grados)')
                legend(lambda_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
                namefig5 = strcat(folderPathFigs, 'props.png') ;
                print(fig5, namefig5,'-dpng') ;


                fig6 = figure(6) ;
                subplot(2,1,1)
                hold on
                plot(r_R, Pn(:, h, i))
                xlabel('r/R')
                ylabel('Fuerza Normal (N)')
                subplot(2,1,2)
                hold on
                plot(r_R, Pt(:, h, i))
                xlabel('r/R')
                ylabel('Fuerza Tangencial (N)')
                legend(lambda_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
                namefig6 = strcat(folderPathFigs, 'forces.png') ;
                print(fig6, namefig6,'-dpng') ;
            end
        end
    end
end

fig7 = figure(7) ;
labelTitle = 'Grafico CpMax' ;
hold on; grid on
for i=1:length(lambda_interp)
    Cp_max(i) = max(Cp_interp(:,i));
end
[maxCp, indexMaxCp] = max(Cp_max);
max_lambda          = lambda_interp(indexMaxCp) ;
plot(lambda_interp, Cp_max, '--r', 'LineWidth',2.0)
hold on
plot(max_lambda, maxCp, 'ok', 'LineWidth', 5.0)
labx = xlabel('\lambda'); laby = ylabel('C_{p,max}');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig7 = strcat(folderPathFigs, 'CpMax.png') ;
print(fig7, namefig7,'-dpng') ;

Coef_labels = {'Lift Coef', 'Drag Coef', 'Momentum Coef'} ;
fig8 = figure(8);
labelTitle = 'Grafico Coef DU21 A17' ;
hold on; grid on
plot(aoa(:,3), clift(:,3), 'r', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,3), cdrag(:,3), 'b', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,3), cmomt(:,3), 'k', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig8 = strcat(folderPathFigs, 'coefDU21_A17.png') ;
print(fig8, namefig8,'-dpng') ;

fig9 = figure(9);
labelTitle = 'Grafico Coef DU25 A17' ;
hold on; grid on
plot(aoa(:,4), clift(:,4), 'r', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,4), cdrag(:,4), 'b', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,4), cmomt(:,4), 'k', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig9 = strcat(folderPathFigs, 'coefDU25_A17.png') ;
print(fig9, namefig9,'-dpng') ;

fig10 = figure(10);
labelTitle = 'Grafico Coef DU30 A17' ;
hold on; grid on
plot(aoa(:,5), clift(:,5), 'r', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,5), cdrag(:,5), 'b', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,5), cmomt(:,5), 'k', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig10 = strcat(folderPathFigs, 'coefDU30_A17.png') ;
print(fig10, namefig10,'-dpng') ;

fig11 = figure(11);
labelTitle = 'Grafico Coef DU35 A17' ;
hold on; grid on
plot(aoa(:,6), clift(:,6), 'r', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,6), cdrag(:,6), 'b', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,6), cmomt(:,6), 'k', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig11 = strcat(folderPathFigs, 'coefDU35_A17.png') ;
print(fig11, namefig11,'-dpng') ;

fig12 = figure(12);
labelTitle = 'Grafico Coef DU40 A17' ;
hold on; grid on
plot(aoa(:,7), clift(:,7), 'r', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,7), cdrag(:,7), 'b', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,7), cmomt(:,7), 'k', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig12 = strcat(folderPathFigs, 'coefD40_A17.png') ;
print(fig12, namefig12,'-dpng') ;

fig13 = figure(13);
labelTitle = 'Grafico Coef NACA64 A17' ;
hold on; grid on
plot(aoa(:,8), clift(:,8), 'r-', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,8), cdrag(:,8), 'b-', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
hold on
plot(aoa(:,8), cmomt(:,8), 'k-', 'LineWidth',2.0); % Ajusta el número de niveles según tus necesidades
labx = xlabel('AoA'); laby = ylabel('Coef');
legend(Coef_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig13 = strcat(folderPathFigs, 'coefNACA64_A17.png') ;
print(fig13, namefig13,'-dpng') ;

disp('Proceso terminado')