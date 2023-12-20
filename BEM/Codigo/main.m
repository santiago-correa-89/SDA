%% Método BEM
clear; 
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

addpath( genpath( [ pwd '/uBEM'] ) );

path          = 'Data/NREL_5MW' ;
file          = '/NRELOffshrBsline5MW_bladeData.txt';

% Abre el archivo en modo de lectura ('r' para lectura)
fid = fopen(fullfile(path, file), 'r');

% Comprueba si el archivo se abrió correctamente
if fid == -1
    error('No se pudo abrir el archivo.');
end

blade = readtable(fullfile(path, file));


%% Datos del problema
% Rhub    = 2.92/2; % Diámetro del Hub
% Rrot    = 61/2  ; % Diámetro del rotor

Rhub    = 4/2        ; % Diámetro del Hub
Rrot    = 116.2674/2 ; % Diámetro del rotor

r       = blade(:,1) ;          % Vector de Radios de la pala
c       = blade(:,2) ;          % Vector de Cuerda de las secciones
theta_t = deg2rad(blade(:,3)) ; % Vector de Ángulo de Twist de la las secciones en radianes
t       = blade(:,4) ;          % Vector de Espesor de las secciones

theta_p = deg2rad(-5:1:5)' ; % Vector de Ángulo de Twist de la las secciones en radianes

V0     = 10;                % m/s
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
a        = zeros(length(r), length(lambda), length(theta_p));
aprime   = zeros(length(r), length(lambda), length(theta_p));
phi      = zeros(length(r), length(lambda), length(theta_p));
alpha    = zeros(length(r), length(lambda), length(theta_p));
Vrel     = zeros(length(r), length(lambda), length(theta_p));
lambda_r = zeros(length(r), length(lambda), length(theta_p));
Cl       = zeros(length(r), length(lambda), length(theta_p));
Cd       = zeros(length(r), length(lambda), length(theta_p));
Pn       = zeros(length(r), length(lambda), length(theta_p));
Pt       = zeros(length(r), length(lambda), length(theta_p));
kiter    = zeros(length(r), length(lambda), length(theta_p));
err      = zeros(length(r), length(lambda), length(theta_p));
M        = zeros(length(lambda), length(theta_p));
Cp       = zeros(length(lambda), length(theta_p));
Ct       = zeros(length(lambda), length(theta_p));

%% Resolución de método BEM
for g = 1:length(theta_p)
    for h = 1:length(lambda)
        clear Ke;
        for i = 1:length(r)
    
            r_i            = r(i);            % Radio de pala
            c_i            = c(i);            % Cuerda de la sección
            theta_t_i      = theta_t(i);      % Angulo de Pitch
            t_i            = t(i);             % Espesor de la sección
            lambda_r(i, h, g) = lambda(h)*r_i/Rrot;    % Cociente Velocidad del aire y velocidad angular en la sección

            [ a(i, h, g) , aprime(i, h, g) , phi(i, h, g) , alpha(i, h, g), Cl(i, h, g), Cd(i, h, g), err(i, h, g) , kiter(i, h, g) ] = BEM_method(B, c_i, r_i, t_i, lambda(h), theta_t_i, theta_p(g), Rrot, Rhub, Tol, 1, 1, factorCorreccion, path );
   
        end

        Vrel(:, h, g)  = sqrt( (V0*(1-a(:, h, g))).^2 + ((lambda_r(:,h , g)*V0.*(1+aprime(:, h, g))).^2 ) );
        
        Ke             = (1/2)*rho*(Vrel(:, h, g).^2).*c(:);
        Pn(:, h, g)    = Ke.*(Cl(:, h, g).*cos(phi(:, h, g)) + Cd(:, h, g).*sin(phi(:, h, g)));
        Pt(:, h, g)    = Ke.*(Cl(:, h, g).*sin(phi(:, h, g)) - Cd(:, h, g).*cos(phi(:, h, g)));
    end
end

%% Variables of BEM to estimate performance of HAWT

for i = 1:length(theta_p)
    for j= 1:length(lambda)
        [ M(j,i) , Cp(j,i) ] = powerFactor(r, Pt(:, j, i), omega(j), rho, V0, Rrot, B);
        [ Ct(j,i) ]          = thrustFactor(r, Pn(:, j, i), rho, V0, Rrot, B);
    end
end

% Crear una cuadrícula más densa para la interpolación
lambda_interp = linspace(min(lambda), max(lambda), length(lambda)*10); % Ajusta el número de puntos según tus necesidades
theta_p_interp = linspace(min(rad2deg(theta_p)), max(rad2deg(theta_p)), length(theta_p)*10); % Ajusta el número de puntos según tus necesidades
[lambda_interp_grid, theta_p_interp_grid] = meshgrid(lambda_interp, theta_p_interp);

% Interpolar los valores de Cp en la cuadrícula densa
Cp_interp = griddata(rad2deg(theta_p), lambda , Cp, theta_p_interp_grid, lambda_interp_grid, 'cubic'); % Puedes ajustar el método de interpolación
Ct_interp = griddata(rad2deg(theta_p), lambda , Ct, theta_p_interp_grid, lambda_interp_grid, 'cubic'); % Puedes ajustar el método de interpolación

% Crear un gráfico de superficie suavizado con contourf
figure(1)
contourf(lambda_interp_grid, theta_p_interp_grid, Cp_interp, 20, '--', "ShowText",true,"LabelFormat", "%0.2f"); % Ajusta el número de niveles según tus necesidades
xlabel('\lambda');
ylabel('\theta_p (grados)');
title('Grafico de superficie Cp');

% Crear un gráfico de superficie suavizado con contourf
figure(2)
contourf(lambda_interp_grid, theta_p_interp_grid, Ct_interp, 20, '--', "ShowText",true,"LabelFormat", "%0.2f"); % Ajusta el número de niveles según tus necesidades
xlabel('\lambda');
ylabel('\theta_p (grados)');
title('Grafico de superficie Ct');

theta_p_labels = {'\theta_p = -5', '\theta_p = -4', '\theta_p = -3', '\theta_p = -2', '\theta_p = -1', '\theta_p = 0', '\theta_p = 1', '\theta_p = 2', '\theta_p = 3', '\theta_p = 4', '\theta_p = 5'};
lambda_labels  = {'\lambda = 7', '\lambda = 9', '\lambda = 11'};


for l =1:length(theta_p)
    figure(3)
    aux = length(Cp_interp)/length(theta_p);
    plot(lambda_interp, Cp_interp((l-1)*aux+1,:))
    hold on
    xlabel('\lambda');
    ylabel('Cp');
    title('Grafico Cp-\lambda');
end

legend(theta_p_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades

for l =1:length(theta_p)
    figure(4)
    aux = length(Ct_interp)/length(theta_p);
    plot(lambda_interp, Ct_interp((l-1)*aux+1,:))
    hold on
    xlabel('\lambda');
    ylabel('Cp');
    title('Grafico Ct-\lambda');
end

legend(theta_p_labels, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades

for i = 1:length(theta_p)
    for h = 1:length(lambda)
        r_R = r/Rrot;
        if theta_p(i) == 0
            if lambda(h) == 7 || lambda(h) == 9 || lambda(h) == 11
                
                figure(5)
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


                figure(6)
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
            end
        end
    end
end

figure(7)
for i=1:length(lambda_interp)
    Cp_max(i) = max(Cp_interp(:,i));
end
plot(lambda_interp, Cp_max, '--r', 'LineWidth',2.0)
ylabel('C_{p, max}')
xlabel('\lambda')
% figure(7)
% for i = 1:length(theta_p)
%     for h = 1:length(lambda)
%         for j = 1:length(r)
%             if theta_p(i) == 0
%                 if r(j) == 31
%                     if lambda(h) == 7
%                         subplot(3,1,1)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 7')
%                     elseif lambda(h) == 9
%                         subplot(3,1,2)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 9')
%                     elseif lambda(h) == 11
%                         subplot(3,1,3)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 11')
%                     end
%                 end
%             end
%         end
%     end
% end

% figure(8)
% for i = 1:length(theta_p)
%     for h = 1:length(lambda)
%         for j = 1:length(r)
%             if theta_p(i) == 0
%                 if r(j) == 54.3
%                     if lambda(h) == 7
%                         subplot(3,1,1)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 7')
%                     elseif lambda(h) == 9
%                         subplot(3,1,2)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 9')
%                     elseif lambda(h) == 11
%                         subplot(3,1,3)
%                         aux = [0, 0];
%                         Vax = [0, (1 - a(j, h, i))*V0];
%                         Vta = [(1 + aprime(j, h, i))*V0*lambda_r(j, h, i), (1 - a(j, h, i))*V0];
%                         plot([aux(1) Vax(1)],[aux(2) Vax(2)],'--^r');
%                         hold on
%                         plot([aux(1) Vta(1)],[Vta(2) Vta(2)],'-->b');
%                         hold on
%                         plot([aux(1) Vta(1)],[aux(2) Vax(2)],'--k');
%                         xlabel('V (m/s)')
%                         ylabel('V (m/s)')
%                         legend('\lambda = 11')
%                     end
%                 end
%             end
%         end
%     end
% end

disp('Proceso terminado')