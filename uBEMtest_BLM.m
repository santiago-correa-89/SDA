clear all ;
close all ;
clc ;
addpath( genpath( [ pwd '/Control'] ) );
addpath( genpath( [ pwd '/Turbsim'] ) );
addpath( genpath( [ pwd '/uBEM']    ) );

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
I      =  nGen*nGen*Igen + Irot + 3*Iblade ; % Total inertia

% Flags
DWM     = false; % Turn on Dynamic Wake Model
DSM     = false; % Turn on Dynamic Stall Model
OyeDSM  = false; % Turn on Oye Dynamic stall model
BLDSM   = false ; % Turn on B-L Dynamic stall model

% Wind Flags
Uniform = false ; % Turn on uniform wind profile with X direction
Share   = true ; % Turn on Share wind profile
Turbsim = false  ; % Read Turbsim file to read wind velocity in the rotor plane

% Wind initialization 
if Turbsim
    FileName = 'ejemploParte2.bts' ;
    [Vwind, twrV0, z, y, yTwr, nz, ny, dy, dz, dt, yHub, y1, meanVhub] = readfile_BTS( FileName ) ;
elseif Share
    Vhub   = 18        ; % Wind velocity at hub height
    n      = 1/7       ; % Wind share potential coef
elseif Uniform
    Vhub   = 15        ;
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
if ~Turbsim
    delta_t    = 0.05  ;
    maxtime    = 150    ;   
    timeVector = (0:delta_t:maxtime) ;
    tIter      = length(timeVector)  ;
elseif Turbsim
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
[aoast, liftCoef, dragCoef, momCoef ] = readairfoildata_v3( polars, path  ) ;
% Data extracted from NREL 5MW documentation
AoA0data  = [  0,    0,    -4.2,    -3.2,    -2,2,    -1.2,    -3.2,  -4.432] ;
AoA1data  = [  0,    0,     8.0,     8.5,       9,    11.5,       9,       9] ;
AoA2data  = [  0,    0,    -8.0,    -8.5,      -9,   -11.5,      -9,      -9] ;
CnAoAdata = [  0,    0,  6.2047,  6.4462,  7.3326,  7.1838,  7.4888,  6.0031] ;
Cd0data   = [0.5, 0.35,   0.006,   0.006,   0.008,   0.012,    0.03,  0.0065] ;
Cn1data   = [  0,    0,  1.4144,  1.4336,   1.449,  1.6717,  1.3519,  1.4073] ;
Cn2data   = [  0,    0, -0.5324, -0.6873, -0.6138, -0.3075, -0.3226, -0.7945] ;

%initialization variables
if Turbsim
    Vx              = zeros(yGridLen, zGridLen,  tIter) ;
    Vy              = zeros(yGridLen, zGridLen,  tIter) ;
    Vz              = zeros(yGridLen, zGridLen,  tIter) ;
end
Vxhub           = zeros( tIter ) ;
Vyhub           = zeros( tIter ) ;
Vzhub           = zeros( tIter ) ;
Vxblade         = zeros( tIter ) ;
Vyblade         = zeros( tIter ) ;
Vzblade         = zeros( tIter ) ;
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
AoA             = zeros(ni, nblade,  tIter) ;
flowangle_deg   = zeros(ni, nblade,  tIter) ;
% BL model init + factors
thetaPitch      = zeros(tIter, nblade)      ;
if BLDSM
    T               = [ 1.5, 5.0, 6.0, 11.0 ]   ;    %  LBM model time constant Tp, Tf0, Tv0, Tvl ( ref Pereira 2011 )
    [ BLcoefs ]     = initDSMbeddoesLeishman( ni, nblade, tIter ) ;
end

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

[ thetaPitch(1,:), lastTimePC(1,1), ~, integError(1,1) ] = initPitchControl( nGen, omega(1, 1), ...
                                                            thetaPitch(1,:), integError(1,1), ...
                                                            lastTimePC(1,1), timeVector(1,1) );

%Initialization ends
%Main loop starts in time domain

for nt = 2:length(timeVector)-1
    timeVector(nt)
    if Turbsim
        % Generate velocity componente in turbsim system for nt
        % time step
        Vx(:,:, nt) = squeeze( Vwind(nt,1,:,:) ) ; Vy(:,:, nt) = squeeze( Vwind(nt,2,:,:) ) ; Vz(:,:, nt) = squeeze( Vwind(nt,3,:,:) );
        % Hub velocity
        yHub = round(yGridLen/2);  zHub = round(zGridLen/2);
        Vxhub(nt) = Vwind(nt, 1, yHub, zHub) ;
        Vyhub(nt) = Vwind(nt, 2, yHub, zHub) ;
        Vzhub(nt) = Vwind(nt, 3, yHub, zHub) ;
    end
    
    if timeVector(nt) >= 10
        DWM   = true ;  % Turn on Dynamic Wake Model
        DSM   = false ;  % Turn on Dynamic Stall Model
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
        
        for j=1:ni

            % Update global position of the blade node rb(i)
            ID            = sectionID(j)      ;                  % Blade section ID
            c             = chord(j)          ;                  % Blade sectio chord
            rblade        = [0, radius(j), 0]';                  % Blade section position vector in local system
            rGb           = a14'*rblade + a12'*Xhub + Ytower ;   % Blade section position vector in global system
            rGx(j, i, nt) = rGb(1);
            rGy(j, i, nt) = rGb(2);
            rGz(j, i, nt) = rGb(3);
            rRot          = a34'*rblade ;     % Blade section position vector in system 3
    
            %Estimate relative wind speed (velocity triangle)            
            if Share
                V0g   = shareWind(t, rGy(j, i, nt), Vhub, Hhub(2), n); % wind velocity in global coordinates share wind
            elseif Uniform
                V0g   = uniWind( t, rGy(j, i, nt), Vhub ) ;            % wind velocity in global coordinates uniform wind
            elseif Turbsim
                % Interpolate coordinates velocity for actual section location
                V0x(j, i, nt) = interp2( Ygrid, Zgrid, Vx(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0y(j, i, nt) = interp2( Ygrid, Zgrid, Vy(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0z(j, i, nt) = interp2( Ygrid, Zgrid, Vz(:,:, nt), rGy(j, i, nt) , rGz(j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0g           = [ V0x(j, i, nt); V0y(j, i, nt); V0z(j, i, nt) ] ;
                Egl           = [ 1 0 0 ; 0 0 -1; 0 1 0];
                V0g           = Egl*V0g ;
                if j == ni-1
                    if i == 1
                        Vxblade(nt) = V0x(j, i, nt);
                        Vyblade(nt) = V0y(j, i, nt);
                        Vzblade(nt) = V0z(j, i, nt);
                    end
                end
            end

            V0l   = a14*V0g ;                                    % convert wind velocity in global coordinate to local blade coordinates            
            L2    = [1 0 0; 0 0 0; 0 0 1];
            V0l   = L2*V0l ;

            W     = [wx(j, i, nt-1) 0 wz(j, i, nt-1)]' ;
            vRot  = [0 0 omega(nt,1)*rRot(2)]';
            vRel  = V0l + W - vRot;  

            [flowangle, vCoordSys ] = computeFlowAngle( vRel ) ;

            vrelNorm = norm(vRel) ;  
            flowangle_deg(j, i, nt)   = atan( vRel(3)/vRel(2)) ; %rad2deg(flowangle)                ; % Convert to degress to interpolate in airfoils coefs    
            thetaTwist                = ( thetaPitch(nt, i) + twist(j) )  ; % Angle in DEGRESS
            AoA(j, i, nt)             = flowangle_deg(j, i, nt) - thetaTwist        ;
            
            %interpolate to find lift and drag (static)      
            %interpolate the values of lift and drag to the different thicknesses                
                       
            clstat    =  interp1( aoast(:,ID), liftCoef(:,ID),  AoA(j, i, nt), 'linear', 'extrap' );
            cdstat    =  interp1( aoast(:,ID), dragCoef(:,ID),  AoA(j, i, nt), 'linear', 'extrap' );
            
            % Apply Stall Delay Model 
            % Apply Dynamic Stall Model
            if DSM
                if OyeDSM
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, ...
                                                                        clfs, delta_t, vrelNorm, chord(j));
                    cdrag(j, i, nt)   = cdstat;
                
                elseif BLDSM
                    if ID == 1 && ID == 2  %% no applying model to cylinder section
                        clift(j, i, nt)   = clstat;
                        cdrag(j, i, nt)   = cdstat;
                    else
                        boolLeishmann = true ;
                        AoA0 = AoA0data(ID);
                        AoA1 = AoA1data(ID);
                        AoA2 = AoA2data(ID);
                        CnAo = CnAoAdata(ID);
                        Cd0  = Cd0data(ID);
                        Cn1  = Cn1data(ID);
                        Cn2  = Cn2data(ID);
        
                        [ BLcoefs(:, j, i, nt) ] = DSMbeddoessLeishman(  AoA(j, i, nt), vrelNorm, delta_t, c, ...
                                                          BLcoefs(:, j, i, nt-1), AoA0,  Cd0, CnAo, ...
                                                          clstat, cdstat, Cn1, Cn2, AoA1, AoA2, T , ...
                                                          boolLeishmann );
                        
                        clift(j, i, nt)  = BLcoefs(1, j, i, nt) ; 
                        cdrag(j, i, nt)  = BLcoefs(2, j, i, nt) ;

                    end
                else
                    clift(j, i, nt)   = clstat;
                    cdrag(j, i, nt)   = cdstat;
                end
            else
                clift(j, i, nt)   = clstat;
                cdrag(j, i, nt)   = cdstat;
            end
           
            %calculate aerodynamic loads
            lift = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*clift(j, i, nt); % Lift force
            drag = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*cdrag(j, i, nt); % Drag force
            
            BEMforces = vCoordSys * [lift, 0 , drag]'; % Tangencial Force local blade system
            
            px(j, i, nt) = BEMforces(1) ;
            pz(j, i, nt) = BEMforces(3) ;

            % For time step update induced wind (relaxation method)
            % Apply glauert heavy induced factor correction
            V0 = norm(V0l) ;
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
            
            if DWM
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
            moment(nt,nblade) = powerFactor( radius(:) , pz(:,nblade,nt) );  
            thrust(nt,nblade) = thrustFactor( radius(:), px(:,nblade,nt) ); 
    end

    aeroPower(nt, 1)   = ( moment(nt, 1) + moment(nt, 2) + moment(nt, 3) )*omega(nt,1); 
    powerOutput(nt, 1) = ( moment(nt, 1) + moment(nt, 2) + moment(nt, 3) )*omega(nt,1)*genEff; %omega(nt);
    totalThrust(nt, 1) = (thrust(nt, 1) + thrust(nt, 2) + thrust(nt, 3))/1000;
    rotTrq(nt, 1)      = moment(nt, 1) + moment(nt, 2) + moment(nt, 3);
    
    %Turn on control after dynamic models are working
    [ trqGen(nt, 1), lastTimeGC(nt, 1), genOmegaF(nt, 1) ] = genControl( nGen, omega(nt, 1), thetaPitch(nt, 1), ...
                               trqGen(nt-1, 1), genOmegaF(nt-1, 1), lastTimeGC(nt-1, 1), timeVector(nt) ) ;
   
    if nt <= length(timeVector) - 1 && timeVector(nt) > 10
        [ thetaPitch(nt+1, :), lastTimePC(nt+1, 1), speedError(nt+1,1), integError(nt+1, 1) ] = pitchControl( nGen, omega(nt, 1), thetaPitch(nt, :), integError(nt, 1), lastTimePC(nt, 1), timeVector(nt) ) ;
        omega(nt+1, 1)  =  omega(nt, 1) + ( ( rotTrq(nt, 1) - (trqGen(nt, 1)*nGen) )/I ) * delta_t ;
    end
end  %end iteration(timesim)

%
studyCase = 'BLM_share_noDSM';
startPlot = 1; 
tStep     = delta_t ;
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ; 
folderPathFigs = './figs/uBEM/Outputs/fullBLtest' ;
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

if Turbsim
    fig6 = figure(6);
    hold on; grid on
    labelTitle = 'Velocidad de viento en hub';
    plot(timeVector(startPlot:end-1), Vxhub(startPlot:end-1, 1), 'r--', 'LineWidth', lw,'markersize',ms )
    hold on
    plot(timeVector(startPlot:end-1), Vyhub(startPlot:end-1, 1), 'b--', 'LineWidth', lw,'markersize',ms )
    hold on
    plot(timeVector(startPlot:end-1), Vzhub(startPlot:end-1, 1), 'k--', 'LineWidth', lw,'markersize',ms )
    labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad (m/s) ');
    set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
    set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
    title(labelTitle)
    namefig6 = strcat(folderPathFigs, ['/hubVelocity=', studyCase,'.jpg'] ) ;
    print(fig6, namefig6,'-dpng') ;
    
    fig7 = figure(7);
    hold on; grid on
    labelTitle = 'Velocidad de viento en hub en punta de pala';
    plot(timeVector(startPlot:end-1), Vxblade(startPlot:end-1, 1), 'r--', 'LineWidth', lw,'markersize',ms )
    hold on
    plot(timeVector(startPlot:end-1), Vyblade(startPlot:end-1, 1), 'b--', 'LineWidth', lw,'markersize',ms )
    hold on
    plot(timeVector(startPlot:end-1), Vzblade(startPlot:end-1, 1), 'k--', 'LineWidth', lw,'markersize',ms )
    labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad (m/s) ');
    set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
    set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
    title(labelTitle)
    namefig7 = strcat(folderPathFigs, ['/bladeVelocity=', studyCase,'.jpg'] ) ;
    print(fig7, namefig7,'-dpng') ;
end

fig8 = figure(8);
hold on; grid on
labelTitle = 'Fuerza de empuje sobre el rotor';
plot(timeVector(startPlot:end-1), totalThrust(startPlot:end-1, 1)', 'r-o', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' Fuerza (N) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig8 = strcat(folderPathFigs, ['/totalThrust=', studyCase,'.jpg'] ) ;
print(fig8, namefig8,'-dpng') ;

cl     = zeros(length(timeVector(1:end-1)), 2);
cd     = zeros(length(timeVector(1:end-1)), 2);
AoAout = zeros(length(timeVector(1:end-1)), 2);
for i = 1:length(timeVector) - 1
    AoAout(i,1:2) = [ AoA(7, 1, i)  , AoA(15, 1, i) ];
    cl(i,1:2)     = [ clift(7, 1, i), clift(15,1,i) ];
    cd(i,1:2)     = [ cdrag(7, 1, i), cdrag(15,1,i) ];   
end

fig9 = figure(9);
hold on; grid on
labelTitle = 'Lift coef';
plot(timeVector(1000:end-1), cl(1000:end,1)', 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1000:end-1), cl(1000:end,2)', 'b', 'LineWidth', lw,'markersize',ms )
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
labx=xlabel( 'Time (s)' ); laby=ylabel(' C_l ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig9 = strcat(folderPathFigs, ['/liftCoef=', studyCase,'.jpg'] ) ;
print(fig9, namefig9,'-dpng') ;

fig10 = figure(10);
hold on; grid on
labelTitle = 'Lift coef';
plot(timeVector(1000:end-1), cd(1000:end,1), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1000:end-1), cd(1000:end,2), 'b ', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' C_d ');
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig10 = strcat(folderPathFigs, ['/liftCoef=', studyCase,'.jpg'] ) ;
print(fig10, namefig10,'-dpng') ;

fig11 = figure(11);
hold on; grid on
labelTitle = 'C_l - AoA';
plot(AoAout(1000:end, 1), cl(1000:end,1), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoAout(1000:end, 2), cl(1000:end,2), 'b ', 'LineWidth', lw,'markersize',ms )
labx=xlabel('AoA (deg)'); laby=ylabel('C_l');
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig11 = strcat(folderPathFigs, ['/liftCoef_AoA=', studyCase,'.jpg'] ) ;
print(fig11, namefig11,'-dpng') ;

fig12 = figure(12);
hold on; grid on
labelTitle = 'AoA';
plot(timeVector(1:end-1), AoAout(:,1), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1:end-1), AoAout(:,2), 'b', 'LineWidth', lw,'markersize',ms )
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
labx=xlabel('Time (s)'); laby=ylabel('AoA (deg)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig12 = strcat(folderPathFigs, ['/AoA=', studyCase,'.jpg'] ) ;
print(fig12, namefig12,'-dpng') ;

clear fig1; clear fig2; clear fig3; clear fig4; clear fig5; clear fig6; clear fig7; clear fig8; clear fig9; clear fig10; clear fig11;