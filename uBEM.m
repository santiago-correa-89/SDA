clear all ;
close all ;
clc ;
addpath( genpath( [ pwd '/Control'] ) );
addpath( genpath( [ pwd '/Turbsim'] ) );
addpath( genpath( [ pwd '/uBEM'] ) );

% General inputs
Dh     =  5          ; % Distance between center of the hub and junction point between tower and nacell (m)
R      =  126/2      ; % Blade Radius (m)
Ht     =  90         ; % Tower Heigth (m)
Rhub   =  3/2        ; % Hub Radius (m)
rho    =  1.225      ; % Flux density
nblade =  3          ; % Number of blades
nGen   =  97.0       ; % Generator ratio speed
genEff =  94.4/100   ; % Generator Efficiency
Igen   =  534.116    ; % Generator inertia kg m^2
Irot   =  115926     ; % Rotor inertia
Iblade =  11776047   ; % Blade inertia in the root 
I      = nGen*nGen*Igen + Irot + 3*Iblade ; % Total inertia

% Flags
DWM   = false; % Turn on Dynamic Wake Model
DSM   = false; % Turn on Dynamic Stall Model
SDM   = false; % Turn on Stall Delay Model
% Wind Flags
Uniform = false ; % Turn on uniform wind profile with X direction
Share   = false ; % Turn on Share wind profile
Turbsim = true  ; % Read Turbsim file to read wind velocity in the rotor plane

% Wind initialization 
if Turbsim == 1
    FileName = 'ejemploParte2.bts' ;
    [Vwind, twrV0, z, y, yTwr, nz, ny, dy, dz, dt, yHub, y1, meanVhub] = readfile_BTS( FileName ) ;
else
    Vhub   = 20      ; % Wind velocity at hub height
    n      = 1/7       ; % Wind share potential coef
end

% Blade geometry and mass distribution
path       = '/Data/NREL_5MW' ;
file       = '/Polars/NRELOffshrBsline5MWbladeDataV1.txt';
bladeData   = readtable(fullfile(path, file));
radius      = bladeData.RNodes    ;
twist       = bladeData.AeroTwst  ; % Angle in DEGRESS
chord       = bladeData.Chord     ;
sectionID   = bladeData.NFoil     ;
ni          = max(size(radius))   ;

Ytower      = [0   , Ht , 0 ]' ;
Xhub        = [-Dh , 0   , 0 ]' ;

thetaTilt   = 5*pi/180 ; % Angle must be in RADIANS!
thetaYaw    = 0*pi/180 ; % Angle must be in RADIANS!
thetaCone   = 2.5*pi/180 ; % Angle must be in RADIANS!

%initialization starts
if Turbsim == 0
    delta_t    = 0.05  ;
    timeVector = (0:delta_t:maxtime) ;
    tIter      = length(timeVector)  ;
elseif Turbsim == 1
    delta_t = dt;
    [tLen, vLen, zGridLen, yGridLen] = size(Vwind) ;
    maxTime    = tLen*delta_t;
    timeVector = (0:delta_t:maxTime);
    tIter      = length(timeVector)  ;
    % Generate grid vector 
    % Creat a grid with the Y and Z coordinates of the rotor
    % plane output of the TurbSim
    [ Ygrid, Zgrid ] = meshgrid(y , z);
end

omega      = zeros(tIter, 1) ;     % Initialize angular speed rad/seg
%omega     = 1.2671           ;
%Read pre-determined airfoil data for the various blade elements.
%The airfoil data is for this turbine described by 6 files containing 
%lift, drag and moment coefficients for 6 different thickness to chord ratio.
%The local airfoil data on the blade are afterwards interpolated to the 
%actual thickness. The thickness and the file names for the six airfoils
%are hardcoded below.

file1      = '/Polars/Cylinder1.txt' ;
file2      = '/Polars/Cylinder2.txt' ;
file3      = '/Polars/DU21_A17.txt'  ;
file4      = '/Polars/DU25_A17.txt'  ;
file5      = '/Polars/DU30_A17.txt'  ;
file6      = '/Polars/DU35_A17.txt'  ;
file7      = '/Polars/DU40_A17.txt'  ;
file8      = '/Polars/NACA64_A17.txt';

polars = { file1, file2, file3, file4, file5, file6, file7, file8 } ;

% check airfoil data of the section using the ID 
[aoa, liftCoef, dragCoef, momCoef ] = readairfoildata_v3( polars, path ) ;

%initialization variables
Vx              = zeros(yGridLen, zGridLen,  tIter) ;
Vy              = zeros(yGridLen, zGridLen,  tIter) ;
Vz              = zeros(yGridLen, zGridLen,  tIter) ;
Vhub            = zeros( tIter ) ;
Vyhub           = zeros( tIter ) ;
Vzhub           = zeros( tIter ) ;
rGx             = zeros(ni, nblade,  tIter) ;
rGy             = zeros(ni, nblade,  tIter) ;
rGz             = zeros(ni, nblade,  tIter) ;
V0x             = zeros(ni, nblade,  tIter) ;
V0y             = zeros(ni, nblade,  tIter) ;
V0z             = zeros(ni, nblade,  tIter) ;
wxqs            = zeros(ni, nblade,  tIter) ;
wzqs            = zeros(ni, nblade,  tIter) ;
wxint           = zeros(ni, nblade,  tIter) ;
wzint           = zeros(ni, nblade,  tIter) ;
wx              = zeros(ni, nblade,  tIter) ;
wz              = zeros(ni, nblade,  tIter) ;
px              = zeros(ni, nblade,  tIter) ;
pz              = zeros(ni, nblade,  tIter) ;
fs              = zeros(ni, nblade,  tIter) ;
clift           = zeros(ni, nblade,  tIter) ;
cdrag           = zeros(ni, nblade,  tIter) ;
angle_of_attack = zeros(ni, nblade,  tIter) ;
flowangle_deg   = zeros(ni, nblade,  tIter) ;
thetaPitch      = zeros(tIter, nblade)      ;

% Torque control values
trqGen            = zeros(tIter, 1) ;
lastTimeGC        = zeros(tIter, 1) ;
genOmegaF         = zeros(tIter, 1) ;
% Pitch control values
integError        = zeros(tIter, 1) ;
speedError        = zeros(tIter, 1) ;
lastTimePC        = zeros(tIter, 1) ;

% Outputs
theta_w     = zeros(tIter, nblade) ;
aeroPower   = zeros(tIter, 1) ;
powerOutput = zeros(tIter, 1) ;
totalThrust = zeros(tIter, 1) ;
rotTrq      = zeros(tIter, 1) ;
moment      = zeros(tIter, nblade) ;
thrust      = zeros(tIter, nblade) ;
Pz1         = zeros(tIter, 2) ;
Px1         = zeros(tIter, 2) ;
Wz1         = zeros(tIter, 2) ;
Wx1         = zeros(tIter, 2) ;
AoAtt       = zeros(tIter, 1) ;
CL          = zeros(tIter, 1) ;

theta_w(1,1) = 0*pi/180          ; % Angle must be in RADIANS!
theta_w(1,2) = theta_w(1)+2*pi/3 ; % Angle must be in RADIANS!
theta_w(1,3) = theta_w(1)+4*pi/3 ; % Angle must be in RADIANS!       

%Init pitch control parameters

    [ thetaPitch(1,:), lastTimePC(1,1), ~, integError(1,1) ] = initPitchControl( nGen, omega(1, 1), thetaPitch(1,:), integError(1,1), ...
                                                            lastTimePC(1,1), timeVector(1,1) );

%Initialization ends
%Main loop starts in time domain

for nt = 2:length(timeVector)-1
    timeVector(nt)
    if Turbsim == 1
        % Generate velocity componente in turbsim system for nt
        % time step
        Vx(:,:, nt) = squeeze( Vwind(nt,1,:,:) ) ; Vy(:,:, nt) = squeeze( Vwind(nt,2,:,:) ) ; Vz(:,:, nt) = squeeze( Vwind(nt,3,:,:) );
        % Hub velocity
        yHub = round(yGridLen/2);  zHub = round(zGridLen/2);
        Vhub(nt)  = Vwind(nt, 1, yHub, zHub) ;
        Vyhub(nt) = Vwind(nt, 2, yHub, zHub) ;
        Vzhub(nt) = Vwind(nt, 3, yHub, zHub) ;
    end
    
    if timeVector(nt) >= 10
        DWM   = true ;  % Turn on Dynamic Wake Model
        DSM   = false ;  % Turn on Dynamic Stall Model
        SDM   = false ;  % Turn on Stall Delay Model
    end

    if nt == 2
            [ thetaPitch(nt,:), lastTimePC(nt,1), ~, integError(nt,1) ] = initPitchControl( nGen, omega(nt, 1), thetaPitch(nt-1,:), integError(nt-1,1), ...
                                                            lastTimePC(nt-1,1), timeVector(nt) );
    end
     
    % Updating blade azimuthel positions
    theta_w(nt, 1) = theta_w(nt-1, 1) + omega(nt,1)*delta_t; % omega
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
            ID            = sectionID(j) ; 
            rblade        = [0, radius(j), 0]';                  % Blade section position vector in local system
            rGb           = a14'*rblade + a12'*Xhub + Ytower ;   % Blade section position vector in global system
            rGx(j, i, nt) = rGb(1);
            rGy(j, i, nt) = rGb(2);
            rGz(j, i, nt) = rGb(3);
            rRot          = a34'*rblade ;     % Blade section position vector in system 3
    
            %Estimate relative wind speed (velocity triangle)
            
            if Share == 1
                V0g   = shareWind(t, rGy(j, i, nt), Vhub, Hhub(2), n); % wind velocity in global coordinates share wind
            elseif Uniform == 1
                V0g   = uniWind( t, rGy(j, i, nt), Vhub ) ;            % wind velocity in global coordinates uniform wind
            elseif Turbsim == 1
                % Interpolate coordinates velocity for actual section location
                V0x(j, i, nt) = interp2( Ygrid, Zgrid, Vx(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0y(j, i, nt) = interp2( Ygrid, Zgrid, Vy(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0z(j, i, nt) = interp2( Ygrid, Zgrid, Vz(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0g = [ V0x(j, i, nt); V0y(j, i, nt); V0z(j, i, nt) ] ;
                Egl = [ 1 0 0 ; 0 0 -1; 0 1 0];
                V0g = Egl*V0g ;
            end

            V0l   = a14*V0g ;                                    % convert wind velocity in global coordinate to local blade coordinates            
            L2    = [1 0 0; 0 0 0; 0 0 1];
            V0l   = L2*V0l ;

            W     = [wx(j, i, nt-1) 0 wz(j, i, nt-1)]' ;
            vRot  = [0 0 omega(nt,1)*rRot(2)]';
            vRel  = V0l + W - vRot;  

            [flowangle, vCoordSys ] = computeFlowAngle( vRel ) ;

            % vrelz    = V0l(3) + wz(j, i, nt-1) - omega*rRot(2) ;  %  Tangecial relavite velocity
            % vrelx    = V0l(1) + wx(j, i, nt-1)                        ;  % Normal relative velocity
            % vrel_sqr = vrelx^2 + vrelz^2;
            % flowangle                 = atan( vrelx/-vrelz )              ;
            vrelNorm = norm(vRel) ;  
            flowangle_deg(j, i, nt)   = rad2deg(flowangle)                ; % Convert to degress to interpolate in airfoils coefs    
            thetaTwist                = ( thetaPitch(nt, i) + twist(j) )  ; % Angle in DEGRESS
            angle_of_attack(j, i, nt) = flowangle_deg(j, i, nt) - thetaTwist        ;
            
            %interpolate to find lift and drag (static)      
            %interpolate the values of lift and drag to the different thicknesses                
          
            clstat    =  interp1( aoa(:,ID), liftCoef(:,ID),  angle_of_attack(j, i, nt), 'linear', 'extrap' );
            cdstat    =  interp1( aoa(:,ID), dragCoef(:,ID),  angle_of_attack(j, i, nt), 'linear', 'extrap' );
            %clinvthick(k) =  interp1( aoa(:,k), clinv(:,k),      angle_of_attack(j, i, nt), 'linear', 'extrap' );
            %clfsthick(k)  =  interp1( aoa(:,k), clfullysep(:,k), angle_of_attack(j, i, nt), 'linear', 'extrap' );
            %fsthick(k)    =  interp1( aoa(:,k), fstat(:,k),      angle_of_attack(j, i, nt), 'linear', 'extrap' );
      
            % Interpolate to the actual thickness
            %clstat   = interp1( thick_prof(:), clthick(:)    , thick(j), 'linear', 'extrap' ) ; % Static 2D  Cl Coef
            %cdstat   = interp1( thick_prof(:), cdthick(:)    , thick(j), 'linear', 'extrap' ) ; % Static 2D  Cd Coef
            %clinvisc = interp1( thick_prof(:), clinvthick(:) , thick(j), 'linear', 'extrap' ) ; % Invicid 2D Cl Coef
            %clfs     = interp1( thick_prof(:), clfsthick(:)  , thick(j), 'linear', 'extrap' ) ; % Fully Separated Cl Coef 
            %fst      = interp1( thick_prof(:), fsthick(:)    , thick(j), 'linear', 'extrap' ) ; % Separation function
            
            if SDM == 1
                %% Apply Stall Delay Model 
                clstSDM = stallDelayModel(clinvisc, clstat, thetaTwist, chord(j), radius(j));
                if DSM == 1
                    if j ~= 1
                        fst     = ( 2 * sqrt( clstSDM / clinvisc ) - 1)^2 ;
                        %% Apply Dynamic Stall Model
                    end
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, clfs, delta_t, vrelNorm, chord(j));
                    cdrag(j, i, nt)   = cdstat;
                else
                    clift(j, i, nt)   = clstSDM;
                    cdrag(j, i, nt)   = cdstat;
                end
            elseif SDM == 0
                % Not apply Stall Delay Model
                if DSM == 1
                    %% Apply Dynamic Stall Model
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, clfs, delta_t, vrelNorm, chord(j));
                    cdrag(j, i, nt)   = cdstat;
                else
                    clift(j, i, nt)   = clstat;
                    cdrag(j, i, nt)   = cdstat;
                end
            end
           
            %calculate aerodynamic loads
            lift = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*clift(j, i, nt); % Lift force
            drag = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*cdrag(j, i, nt); % Drag force
            
            BEMforces = vCoordSys * [lift, 0 , drag]'; % Tangencial Force local blade system
            
            px(j, i, nt) = BEMforces(1) ;
            pz(j, i, nt) = BEMforces(3) ;
            %[thetaPitch(nt, i), flowangle_deg(j, i, nt), angle_of_attack(j, i, nt), px(j, i, nt), pz(j, i, nt) ]
            
            if isnan(flowangle_deg) == true
                disp('error')
            end
            
            % For time step update induced wind (relaxation method)
            % Apply glauert heavy induced factor correction
            V0 = norm(V0l) ;
            %V0 = sqrt(V0l(1)^2 + V0l(3)^2);
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
                wxint(j, i, nt-1), wz(j, i, nt-1), wzqs(j, i, nt), wzqs(j, i, nt-1), wzint(j, i, nt-1), vrelNorm, radius(j), R, a, delta_t );
            else
                wx(j, i, nt) = wxqs(j, i, nt);
                wz(j, i, nt) = wzqs(j, i, nt);
            end
      
        end   %end blade section 
    end   %end number of blades

    % Calculate thrust and power   
    for nblade=1:3 
            moment(nt,nblade) = powerFactor( radius(:) , pz(:,nblade,nt) ) ; 
            thrust(nt,nblade) = thrustFactor( radius(:), px(:,nblade,nt) ) ;
    end

    aeroPower(nt, 1)   = ( moment(nt, 1) + moment(nt, 2) + moment(nt, 3) )*omega(nt,1); 
    powerOutput(nt, 1) = ( moment(nt, 1) + moment(nt, 2) + moment(nt, 3) )*omega(nt,1)*genEff; %omega(nt);
    totalThrust(nt, 1) = (thrust(nt, 1) + thrust(nt, 2) + thrust(nt, 3))/1000;
    rotTrq(nt, 1)      = moment(nt, 1) + moment(nt, 2) + moment(nt, 3);
    
    %Turn on control after dynamic models are working
    [ trqGen(nt, 1), lastTimeGC(nt, 1), genOmegaF(nt, 1) ] = genControl( nGen, omega(nt, 1), thetaPitch(nt, 1), ...
                               trqGen(nt-1, 1), genOmegaF(nt-1, 1), lastTimeGC(nt-1, 1), timeVector(nt) ) ;
   
    if nt <= length(timeVector) - 1 
        [ thetaPitch(nt+1, :), lastTimePC(nt+1, 1), speedError(nt+1,1), integError(nt+1, 1) ] = pitchControl( nGen, omega(nt, 1), thetaPitch(nt, :), integError(nt, 1), lastTimePC(nt, 1), timeVector(nt) ) ;
        omega(nt+1, 1)  =  omega(nt, 1) + ( ( rotTrq(nt, 1) - (trqGen(nt, 1)*nGen) )/I ) * delta_t ;
    end
end  %end iteration(timesim)

%
studyCase = 'TurbSimParte1';
startPlot = 1; 
tStep     = 0.05 ;
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ; 
folderPathFigs = './figs/uBEM/Outputs/Turbsim/Parte1' ;
mkdir(folderPathFigs) 

fig1 = figure(1);
labelTitle = ' Rotor power ' ;
hold on; grid on
plot(timeVector(startPlot:end-1), powerOutput(startPlot:end-1,1)', 'r-', 'LineWidth',lw,'markersize',ms );
hold on
plot(timeVector(startPlot:end-1), ones(1,length(timeVector(startPlot:end-1)))*5e6, 'b-', 'LineWidth', lw,'markersize',ms )
labx=xlabel('Time (s)'); laby=ylabel('Power (W)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig1 = strcat(folderPathFigs, [ '/rotorPower_', studyCase,'.jpg' ]) ;
print(fig1, namefig1,'-dpng') ;

fig2 = figure(2);
labelTitle =  ' Aerodynamic Power' ;
hold on; grid on
plot(timeVector(startPlot:end-1), rotTrq(startPlot:end-1)', 'r-', 'LineWidth',lw,'markersize',ms )
hold on
plot(timeVector(startPlot:end-1), aeroPower(startPlot:end-1)', 'b-', 'LineWidth',lw,'markersize',ms )
legend('Rotor Momentum', 'Aerodynamic power', 'location','Best')
labx=xlabel('Time (s)'); laby=ylabel('Momentum (Nm)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig2 = strcat(folderPathFigs, ['/aeroDynamicOutput=', studyCase,'.jpg']) ;
print(fig2, namefig2,'-dpng') ;

fig3 = figure(3);
hold on; grid on
labelTitle = 'Variación velocidad angular rotor';
plot(timeVector(startPlot:end-1), omega(startPlot:end-1, 1)*60/2/pi', 'r-o', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(startPlot:end-1), ones(1,length(timeVector(startPlot:end-1)))*12.1, 'b-', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad angular rotor (rad/seg) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig3 = strcat(folderPathFigs, ['/omegaRotor=', studyCase,'.jpg'] ) ;
print(fig3, namefig3,'-dpng') ;

fig4 = figure(4);
hold on; grid on
labelTitle = 'Variación Ángulo de Pitch Pala 1';
plot(timeVector(startPlot:end-1), thetaPitch(startPlot:end-1, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' \theta_{pitch} (degrees) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig4 = strcat(folderPathFigs, ['/pitchControl_Vhub=', studyCase,'.jpg'] ) ;
print(fig4, namefig4,'-dpng') ;

fig5 = figure(5);
hold on; grid on
labelTitle = 'Error acumulado';
plot(timeVector(1:end-1), integError(1:end-1, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1:end-1), speedError(1:end-1, 1)', 'b-', 'LineWidth', lw,'markersize',ms )
legend('Cumulative Error', 'Error Variation', 'location','Best')
labx=xlabel( 'Time (s)' ); laby=ylabel(' Integral error ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig5 = strcat(folderPathFigs, ['/controlError=', studyCase,'.jpg'] ) ;
print(fig5, namefig5,'-dpng') ;

fig6 = figure(6);
hold on; grid on
labelTitle = 'Velocidad de viento en hub';
plot(timeVector(startPlot:end-1), Vhub(startPlot:end-1, 1), 'r--', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(startPlot:end-1), Vyhub(startPlot:end-1, 1), 'b--', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(startPlot:end-1), Vzhub(startPlot:end-1, 1), 'k--', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad (m/s) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig6 = strcat(folderPathFigs, ['/hubVelocidy=', studyCase,'.jpg'] ) ;
print(fig6, namefig6,'-dpng') ;

fig7 = figure(7);
hold on; grid on
labelTitle = 'Fuerza de empuje sobre el rotor';
plot(timeVector(startPlot:end-1), totalThrust(startPlot:end-1, 1)', 'r-o', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' Fuerza (N) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig7 = strcat(folderPathFigs, ['/totalThrust=', studyCase,'.jpg'] ) ;
print(fig7, namefig7,'-dpng') ;

clear fig1; clear fig2; clear fig3; clear fig4; clear fig5; clear fig6; clear fig7;