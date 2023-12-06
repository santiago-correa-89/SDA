clear all;
close all ;
clc ;

folder = pwd ;
addpath( genpath( folder ))  ;

% General inputs

Dh     =  7.1   ; % Distance between center of the hub and junction point between tower and nacell (m)
R      =  89.17 ; % Blade Radius (m)
Ht     =  119   ; % Tower Heigth (m)
Rhub   =  5.6   ; % Hub Radius (m)
rho    =  1.225 ; % Flux density
nblade =  3     ; % Number of blades

% Wind initialization 
Vhub   = 8.0    ; % Wind velocity at hub height
n      = 1/7    ; % Wind share potential coef

% Flags

DWM   = false; % Turn on Dynamic Wake Model
DSM   = false; % Turn on Dynamic Stall Model
SDM   = false; % Turn on Stall Delay Model
Share = true;  % Turn on Share wind profile
Pitch = false;  % Turn on Pitch variation

% Blade geometry and mass distribution
strucdata   = importdata('bladedat.txt');
radius      = strucdata(:,1);
twist       = strucdata(:,2); % Angle in DEGRESS
chord       = strucdata(:,3);
thick       = strucdata(:,4);
ni          = max(size(radius));

Ytower      = [0   , Ht , 0 ]' ;
Xhub        = [-Dh , 0   , 0 ]' ;

thetaPitch  = 0 ; % Angle must be in RADIANS!
thetaTilt   = 3*pi/180 ; % Angle must be in RADIANS!
thetaYaw    = 4*pi/180 ; % Angle must be in RADIANS!
thetaCone   = 4*pi/180 ; % Angle must be in RADIANS!

%initialization starts
omega      = 0.673; % Radians per second
delta_t    = 0.1;
maxtime    = 200;
timeVector = (0:delta_t:maxtime);
tIter      = length(timeVector);
startPlot  = 500;
tStep      = 1;

%Read pre-determined airfoil data for the various blade elements.
%The airfoil data is for this turbine described by 6 files containing 
%lift, drag and moment coefficients for 6 different thickness to chord ratio.
%The local airfoil data on the blade are afterwards interpolated to the 
%actual thickness. The thickness and the file names for the six airfoils
%are hardcoded below.

thick_prof(1) = 100;
thick_prof(2) = 60;
thick_prof(3) = 48;
thick_prof(4) = 36;
thick_prof(5) = 30.1;
thick_prof(6) = 24.1;

[aoa, cl, cd, cm, fstat, clinv, clfullysep]=...
 readairfoildata_v3('cylinder_ds.txt','FFA-W3-600_ds.txt',...
                    'FFA-W3-480_ds.txt','FFA-W3-360_ds.txt',...
                    'FFA-W3-301_ds.txt','FFA-W3-241_ds.txt');

%initialization variables
rGx             = zeros(ni, nblade,  tIter);
rGy             = zeros(ni, nblade,  tIter);
rGz             = zeros(ni, nblade,  tIter);
wxqs            = zeros(ni, nblade,  tIter);
wzqs            = zeros(ni, nblade,  tIter);
wxint           = zeros(ni, nblade,  tIter);
wzint           = zeros(ni, nblade,  tIter);
wx              = zeros(ni, nblade,  tIter);
wz              = zeros(ni, nblade,  tIter);
px              = zeros(ni, nblade,  tIter);
pz              = zeros(ni, nblade,  tIter);
fs              = zeros(ni, nblade,  tIter);
angle_of_attack = zeros(ni, nblade,  tIter);
clift           = zeros(ni, nblade,  tIter);

% Outputs
theta_w = zeros(tIter, nblade) ;
power   = zeros(tIter, 1);
axForce = zeros(tIter, 1);
totMom  = zeros(tIter, 1);
moment  = zeros(tIter, nblade) ;
thrust  = zeros(tIter, nblade) ;
Pz1     = zeros(tIter, 2);
Px1     = zeros(tIter, 2);
Wz1     = zeros(tIter, 2);
Wx1     = zeros(tIter, 2);
AoA     = zeros(tIter, 1);
CL      = zeros(tIter, 1);

theta_w(1,1) = 0*pi/180          ; % Angle must be in RADIANS!
theta_w(1,2) = theta_w(1)+2*pi/3 ; % Angle must be in RADIANS!
theta_w(1,3) = theta_w(1)+4*pi/3 ; % Angle must be in RADIANS!       

%Initialization ends

%Main loop starts in time domain

for nt = 2:length(timeVector)

    if timeVector(nt) >= 10
        DWM   = true;  % Turn on Dynamic Wake Model
        DSM   = true;  % Turn on Dynamic Stall Model
        SDM   = true;  % Turn on Stall Delay Model
    end

    % Function that change the theta pitch of the system
    if Pitch == 1
        if timeVector(nt) < 100
            thetaPitch  = 0 ; 
        elseif timeVector(nt) >= 100 && timeVector(nt) < 150
            if thetaPitch < 2
                thetaPitch  = thetaPitch + 1*delta_t ; % Angle must be in DEGREES!
            elseif thetaPitch == 2
                thetaPitch = theetaPitch;
            end
        else
            if thetaPitch > 0 
                thetaPitch  = thetaPitch - 1*delta_t;
            elseif thetaPitch == 0
                thetaPitch = thetaPitch;
            end
        end
    end
     
    % Updating blade azimuthel positions
    theta_w(nt, 1) = theta_w(nt-1, 1) + omega*delta_t;
    theta_w(nt, 2) = theta_w(nt, 1)   + 2*pi/3;
    theta_w(nt, 3) = theta_w(nt, 1)   + 4*pi/3;
    % Evluation time step
    t = timeVector(nt); 
    
    % For rigid systems these transformatio matrix system don`t change over time (this can be modify for system with yaw misalignment or tower deflection)
    a1   = yawMatrix( thetaYaw );
    a2   = tiltMatrix( thetaTilt );
    a3   = eye(3,3);
    a12  = (a3*a2)*a1;
    a34  = coneMatrix( thetaCone );

    Hhub = Ytower + a12'*Xhub ; % Hub height in global coordiantes

    %For each element update induced wind and calculate aerodynamic loads
    for i=1:3
        % Time dependecy tranformation matrix
        a23 = azimutalMatrix( theta_w(nt, i) );
        a14 = a34 * a23 * a12; % Transformation matrix from system 1 to 4.
        a41 = a14';
        
        for j=1:ni-1

            % Update global position of the blade node rb(i)

            rblade        = [0, radius(j), 0]';                  % Blade section position vector in local system
            rGb           = a14'*rblade + a12'*Xhub + Ytower ;   % Blade section position vector in global system
            rGx(j, i, nt) = rGb(1);
            rGy(j, i, nt) = rGb(2);
            rGz(j, i, nt) = rGb(3);
            rRot          = a34'*rblade ;     % Blade section position vector in system 3
    
            %Estimate relative wind speed (velocity triangle)
            
            if Share == 1
                V0g   = shareWind(t, rGy(j, i, nt), Vhub, Hhub(2), n); % wind velocity in global coordinates share wind
            else
                V0g   = uniWind( t, rGy(j, i, nt), Vhub ) ;            % wind velocity in global coordinates uniform wind
            end

            V0l   = a14*V0g ;                                    % convert wind velocity in global coordinate to local blade coordinates
   
            vrelz    = V0l(3) + wz(j, i, nt-1) - omega*rRot(2) ;  % Tangecial relavite velocity
            vrelx    = V0l(1) + wx(j, i, nt-1)                 ;  % Normal relative velocity
             
            vrel_sqr = vrelx^2 + vrelz^2;
            
            flowangle                 = atan( -vrelx/vrelz )       ;
            flowangle_deg             = rad2deg(flowangle)         ; % Convert to degress to interpolate in airfoils coefs    
            thetaTwist                = ( thetaPitch + twist(j) )  ; % Angle in DEGRESS
            angle_of_attack(j, i, nt) = flowangle_deg - thetaTwist ;

            %interpolate to find lift and drag (static)
      
            %interpolate the values of lift and drag to the different thicknesses
            for k=1:6                 
          
                clthick(k)    =  interp1( aoa(:,k), cl(:,k),         angle_of_attack(j, i, nt), 'linear', 'extrap' );
                cdthick(k)    =  interp1( aoa(:,k), cd(:,k),         angle_of_attack(j, i, nt), 'linear', 'extrap' );
                clinvthick(k) =  interp1( aoa(:,k), clinv(:,k),      angle_of_attack(j, i, nt), 'linear', 'extrap' );
                clfsthick(k)  =  interp1( aoa(:,k), clfullysep(:,k), angle_of_attack(j, i, nt), 'linear', 'extrap' );
                fsthick(k)    =  interp1( aoa(:,k), fstat(:,k),      angle_of_attack(j, i, nt), 'linear', 'extrap' );
      
            end
      
            % Interpolate to the actual thickness
            clstat   = interp1( thick_prof(:), clthick(:)    , thick(j), 'linear', 'extrap' ) ; % Static 2D  Cl Coef
            cdstat   = interp1( thick_prof(:), cdthick(:)    , thick(j), 'linear', 'extrap' ) ; % Static 2D  Cd Coef
            clinvisc = interp1( thick_prof(:), clinvthick(:) , thick(j), 'linear', 'extrap' ) ; % Invicid 2D Cl Coef
            clfs     = interp1( thick_prof(:), clfsthick(:)  , thick(j), 'linear', 'extrap' ) ; % Fully Separated Cl Coef 
            fst      = interp1( thick_prof(:), fsthick(:)    , thick(j), 'linear', 'extrap' ) ; % Separation function
            
            if SDM == 1
                %% Apply Stall Delay Model 
                clstSDM = stallDelayModel(clinvisc, clstat, thetaTwist, chord(j), radius(j));
                if DSM == 1
                    if j ~= 1
                        fst     = ( 2 * sqrt( clstSDM / clinvisc ) - 1)^2 ;
                        %% Apply Dynamic Stall Model
                    end
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, clfs, delta_t, vrel_sqr, chord(j));
                    cdrag   = cdstat;
                else
                    clift(j, i, nt)   = clstSDM;
                    cdrag   = cdstat;
                end
            elseif SDM == 0
                % Not apply Stall Delay Model
                if DSM == 1
                    %% Apply Dynamic Stall Model
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, clfs, delta_t, vrel_sqr, chord(j));
                    cdrag   = cdstat;
                else
                    clift(j, i, nt)   = clstat;
                    cdrag   = cdstat;
                end
            end
           
            %calculate aerodynamic loads
            lift = 0.5*rho*vrel_sqr*chord(j)*clift(j, i, nt); % Lift force
            drag = 0.5*rho*vrel_sqr*chord(j)*cdrag;           % Drag force
            
            px(j, i, nt) = lift*cos(flowangle) + drag*sin(flowangle); % Normal Force local blade system
            pz(j, i, nt) = lift*sin(flowangle) - drag*cos(flowangle); % Tangencial Force local blade system
            
            % For time step update induced wind (relaxation method)
            % Apply glauert heavy induced factor correction
            V0 = sqrt(V0l(1)^2 + V0l(3)^2);
            a  = -wx(j, i, nt-1)/V0;
            if(a<0.333)
                fglau=1;
            else
                fglau=0.25*(5-3*a);
            end
            
            %% Apply prandtl tip and hub correction factor
            F = prandtlFactor(radius(j), R, Rhub, flowangle);

            %% Update induced wind
            aux1  = sqrt( V0l(3)^2 + ( V0l(1) + fglau*wx(j, i, nt-1) )^2 );
            aux2  = aux1*4*pi*rho*radius(j)*F;
            
            wzqs(j, i, nt)  = -3*lift*sin(flowangle)/aux2;
            wxqs(j, i, nt)  = -3*lift*cos(flowangle)/aux2;
            
            if DWM == 1
                [wx(j, i, nt), wz(j, i, nt), wxint(j, i, nt), wzint(j, i, nt)] = dynamicInflow(wx(j, i, nt-1), wxqs(j, i, nt), wxqs(j, i, nt-1), ...
                wxint(j, i, nt-1), wz(j, i, nt-1), wzqs(j, i, nt), wzqs(j, i, nt-1), wzint(j, i, nt-1), vrel_sqr, radius(j), R, a, delta_t );
            else
                wx(j, i, nt) = wxqs(j, i, nt);
                wz(j, i, nt) = wzqs(j, i, nt);
            end
      
        end   %end blade sections

        % Yawed inflow
%         if thetaYaw ~= 0
%             target_value = 0.7;
%             normalized_values = radius(:) / R;
%             [~, index] = min(abs(normalized_values - target_value));
%             W = [wx(index, i, nt), 0, wz(index, i, nt)]';
%             chi = skewAngle(V0g, W, a34, a12) ; 
%             for j=1:ni-1
%                 [wx(j, i, nt), wz(j, i, nt)] = yawedCorrection(wx(j, i, nt), wz(j, i, nt), chi, theta_w(nt, i), thetaYaw, radius(j), R) ;
%             end
%         end

    end   %end number of blades

% Calculate thrust and power
    
%Initialize thrust and moment

    for nblade=1:3 
            moment(nt,nblade) = powerFactor( radius(:) , pz(:,nblade,nt) ) ; 
            thrust(nt,nblade) = thrustFactor( radius(:), px(:,nblade,nt) ) ;
    end
%     moment = 0;
%     thrust = 0;
%     for nblade=1:3
%         for i=1:ni-1 %
%             moment = moment + 0.5*(radius(i+1)*pz(i+1,nblade,nt) + radius(i)*pz(i,nblade,nt))*...
%                       (radius(i+1)-radius(i));
%             thrust = thrust + 0.5*(px(i+1,nblade,nt)+px(i+1,nblade,nt))*(radius(i+1)-radius(i));
%         end
%     end

%     totMom(nt)  = moment       ;
%     power(nt)   = moment*omega ;
%     axForce(nt) = thrust/1000  ;

power(nt)   = ( moment(nt,1) + moment(nt,2) + moment(nt,3) )*omega;
axForce(nt) = (thrust(nt,1) + thrust(nt,2) + thrust(nt,3))/1000;
totMom(nt)  = moment(nt,1) + moment(nt,2) + moment(nt,3);

end  %end iteration(timesim)

rn = [5, 14] ;
for n = 1:2
    for i = 1:length(timeVector)
        Pz1(i,n)   =  pz(rn(n), 1, i);
        Px1(i,n)   =  pz(rn(n), 1, i);
        Wz1(i,n)   =  wz(rn(n), 1, i);
        Wx1(i,n)   =  wx(rn(n), 1, i);
        AoA(i,n)   =  angle_of_attack(rn(n), 1, i);
        CL(i,n)    =  clift(rn(n), 1, i);
    end
end

%
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
fig1 = figure(1);
hold on; grid on
fig1 = plot(timeVector(startPlot:tStep:end), power(startPlot:tStep:end)', 'r-', 'LineWidth',lw,'markersize',ms )
xlabel('Time (s)')
ylabel('Power (W)')
title('Potencia generada')
saveas(fig1, './Figuras/Power_noDWI_noDSM.jpg');

fig2 = figure(2)
hold on; grid on
plot(timeVector(startPlot:tStep:end), totMom(startPlot:tStep:end)', 'r-', 'LineWidth',lw,'markersize',ms )
xlabel('Time (s)')
ylabel('Momento (Nm)')
title('Variación de Momento Total')
saveas(fig2, './Figuras/Momento_noDWI_noDSM.jpg');

fig3 = figure(3)
hold on; grid on
plot(timeVector(startPlot:tStep:end), Pz1(startPlot:tStep:end, 1)', 'r-', 'LineWidth',lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), Pz1(startPlot:tStep:end, 2)', 'b-', 'LineWidth',lw,'markersize',ms )
xlabel('Time')
ylabel(' Force (kN) ')
title('Fuerza tangencial')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig3, './Figuras/Tangencial_noDWI_noDSM.jpg');
    
fig4 = figure(4)
hold on; grid on
plot(timeVector(startPlot:tStep:end), CL(startPlot:tStep:end, 1)', 'r-', 'LineWidth',lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), CL(startPlot:tStep:end, 2)', 'b-', 'LineWidth',lw,'markersize',ms  )
xlabel('Time')
ylabel(' C_{L} ')
title('Ceoficiente de Sustenentación')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig4, './Figuras/CL_noDWI_noDSM.jpg');

fig5 = figure(5)
hold on; grid on
plot(timeVector(startPlot:tStep:end), AoA(startPlot:tStep:end, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), AoA(startPlot:tStep:end, 2)', 'b-', 'LineWidth', lw,'markersize',ms )
xlabel('Time')
ylabel(' AoA ')
title('Ángulo de Ataque')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig5, './Figuras/AoA_noDWI_noDSM.jpg');

fig6 = figure(6)
hold on; grid on
plot(timeVector(startPlot:tStep:end), Px1(startPlot:tStep:end, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), Px1(startPlot:tStep:end, 2)', 'b-', 'LineWidth', lw,'markersize',ms )
xlabel('Time')
ylabel(' Force (kN) ')
title('Fuerza normal')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig6, './Figuras/Empuje_noDWI_noDSM.jpg');

fig7 = figure(7)
hold on; grid on
plot(AoA(500:1:end, 1), CL(500:1:end, 1), 'r-', 'LineWidth', lw,'markersize',ms )
hold on; grid on
plot(AoA(500:1:end, 2), CL(500:1:end, 2), 'b-', 'LineWidth', lw,'markersize',ms )
xlabel('AoA')
ylabel(' C_{L} ')
title('Sustentación Dinámica')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig7, './Figuras/CL_AoA_noDWI_noDSM.jpg');

fig8 = figure(8)
hold on; grid on
plot(timeVector(startPlot:tStep:end), Wx1(startPlot:tStep:end, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), Wx1(startPlot:tStep:end, 2)', 'b-', 'LineWidth', lw,'markersize',ms )
xlabel('Time (s)')
ylabel(' Velocidad inducida componente X (m/s) ')
title('Variación velocidad inducida axial')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig8, './Figuras/Wx_noDWI_noDSM.jpg');

fig9 = figure(9)
hold on; grid on
plot(timeVector(startPlot:tStep:end), Wz1(startPlot:tStep:end, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
hold on; grid on
plot(timeVector(startPlot:tStep:end), Wz1(startPlot:tStep:end, 2)', 'b-', 'LineWidth', lw,'markersize',ms )
xlabel('Time (s)')
ylabel(' Velocidad inducida componente z (m/s) ')
title('Variación velocidad inducida tangencial')
legend('Radio = 32.3 m', 'Radio = 82.7 m', 'location','Best')
saveas(fig9, './Figuras/Wz_noDWI_noDSM.jpg');


% legend(position, 'Location', 'Best'); % Ajusta la ubicación de la leyenda según tus necesidades