clear all ;
close all ;
clc ;
addpath( genpath( [ pwd '/Control'] ) );
addpath( genpath( [ pwd '/Turbsim'] ) );
addpath( genpath( [ pwd '/uBEM']    ) );

check = 0;

% Flags
DWM       = false;  % Turn on Dynamic Wake Model
DSM       = false; % Turn on Dynamic Stall Model
OyeDSM    = false; % Turn on Oye Dynamic stall model
BLDSM     = false ; % Turn on B-L Dynamic stall model
controlON = true ; % Turn on control pitch and generator

% Wind Flags
Uniform = false ; % Turn on uniform wind profile with X direction
Share   = false ; % Turn on Share wind profile
Turbsim = true  ; % Read Turbsim file to read wind velocity in the rotor plane

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


% Blade geometry and mass distribution
path       = '/Data/NREL_5MW' ;
file       = '/Polars/NRELOffshrBsline5MWbladeDataV1.txt';
bladeData   = readtable(fullfile(path, file));
radius      = bladeData.RNodes    ;
twist       = bladeData.AeroTwst  ; % Angle in DEGRESS
chord       = bladeData.Chord     ;
sectionID   = bladeData.NFoil     ;
ni          = max(size(radius))   ;

Xtower      = [Ht  , 0   , 0   ]' ;
Zhub        = [0   , 0   , -Dh ]' ;

thetaTilt   = 5*pi/180 ; % Angle must be in RADIANS!
thetaYaw    = 0*pi/180 ; % Angle must be in RADIANS!
thetaCone   = 2.5*pi/180 ; % Angle must be in RADIANS!

file  = {'/Polars/Cylinder1.txt', '/Polars/Cylinder2.txt', '/Polars/DU21_A17.txt', '/Polars/DU25_A17.txt',...
         '/Polars/DU30_A17.txt' , '/Polars/DU35_A17.txt' , '/Polars/DU40_A17.txt', '/Polars/NACA64_A17.txt'};

% check airfoil data of the section using the ID 
[aoast, liftCoef, dragCoef, momCoef ] = readairfoildata( path, file  ) ;
% Data extracted from NREL 5MW documentation
AoA0data  = [  0,    0,    -4.2,    -3.2,    -2,2,    -1.2,    -3.2,  -4.432] ;
AoA1data  = [  0,    0,     8.0,     8.5,       9,    11.5,       9,       9] ;
AoA2data  = [  0,    0,    -8.0,    -8.5,      -9,   -11.5,      -9,      -9] ;
CnAoAdata = [  0,    0,  6.2047,  6.4462,  7.3326,  7.1838,  7.4888,  6.0031] ;
Cd0data   = [0.5, 0.35,   0.006,   0.006,   0.008,   0.012,    0.03,  0.0065] ;
Cn1data   = [  0,    0,  1.4144,  1.4336,   1.449,  1.6717,  1.3519,  1.4073] ;
Cn2data   = [  0,    0, -0.5324, -0.6873, -0.6138, -0.3075, -0.3226, -0.7945] ;

% Wind initialization 
if Turbsim
    FileName = 'ejemploParte2.bts' ;
    [Vwind, twrV0, y, x, yTwr, nx, ny, dx, dy, dt, xHub, x1, meanVhub] = readfile_BTS( FileName ) ;
elseif Share
    Vhub   = 15        ; % Wind velocity at hub height
    nPot   = 1/7       ; % Wind share potential coef
elseif Uniform
    Vhub   = 15        ;
end

% time and mapping initialization for different wind conditions
if ~Turbsim
    delta_t    = 0.05  ;
    maxtime    = 150    ;   
    timeVector = (0:delta_t:maxtime) ;
    tIter      = length(timeVector)  ;
elseif Turbsim
    delta_t    = dt;
    [tLen, vLen, yGridLen, xGridLen] = size(Vwind) ;
    maxTime    = tLen*delta_t;
    timeVector = (0:delta_t:maxTime);
    tIter      = length(timeVector)  ;
    % Generate grid vector 
    % Creat a grid with the Y and Z coordinates of the rotor
    % plane output of the TurbSim
    [ Xgrid, Ygrid ] = meshgrid(x , y);
end

% initalize Variables
[ Vblade, rG, V0, thetaWing, omega, WakeQS, WakeInt, Wake, fld, clift, cdrag, AoA, flowAngle, fs ] = initVariables(ni, ...
                                                                                                       nblade, tIter) ;
if Turbsim
    Vhub            = zeros( 3, tIter ) ;
    Vx              = zeros(xGridLen, yGridLen,  tIter) ;
    Vy              = zeros(xGridLen, yGridLen,  tIter) ;
    Vz              = zeros(xGridLen, yGridLen,  tIter) ;
end

if BLDSM
    T               = [ 1.5, 5.0, 6.0, 11.0 ]   ;    %  LBM model time constant Tp, Tf0, Tv0, Tvl ( ref Pereira 2011 )
    [ BLcoefs ]     = initDSMbeddoesLeishman( ni, nblade, tIter ) ;
end

% initialize control Variables
if controlON
    [ genOmega, lastTimeGC, thetaPitch, integError, speedError, lastTimePC ] = initControlVariables(nblade, tIter) ;
    
    [ thetaPitch(1,:), lastTimePC(1,1), ~, integError(1,1) ] = initPitchControl( nGen, omega(1, 1), ...
                                                               thetaPitch(1,:), integError(1,1)   , ...
                                                               lastTimePC(1,1), timeVector(1,1)  ) ;
end

% initialize output variables
[ aeroPower, genPower, blThrust, thrust, genTrq, aeroTrq, rotTrq, CL ] = initOutputVariables(tIter, nblade) ;

omega(1,1)      = 12.1*2*pi/60;
thetaPitch(1,:) = 1*[1, 1, 1 ]';
 
for nt = 2:length(timeVector)-1
    if Turbsim
        % Generate velocity componente in turbsim system for nt
        % time step
        Vx(:,:,nt) = squeeze( Vwind(nt,1,:,:) ) ; Vy(:,:, nt) = squeeze( Vwind(nt,2,:,:) ) ; Vz(:,:, nt) = squeeze( Vwind(nt,3,:,:) );
        % Hub velocity
        xHub = round(xGridLen/2);  yHub = round(yGridLen/2);
        Vhub(:, nt) = Vwind(nt, 1:3, xHub, yHub) ;
    end
    
    if timeVector(nt) >= 20
        DWM   = true ;  % Turn on Dynamic Wake Model
        DSM   = true ;  % Turn on Dynamic Stall Model
    end

    if nt == 2
            [ thetaPitch(nt,:), lastTimePC(nt,1), ~, integError(nt,1) ] = initPitchControl( nGen, omega(nt, 1), ...
                                  thetaPitch(nt-1,:), integError(nt-1,1), lastTimePC(nt-1,1), timeVector(nt) );
    end
     
    % Updating blade azimuthel positions
    thetaWing(nt, 1) = thetaWing(nt-1, 1) + omega(nt,1)*delta_t; % omega
    thetaWing(nt, 2) = thetaWing(nt, 1)   + 2*pi/3;
    thetaWing(nt, 3) = thetaWing(nt, 1)   + 4*pi/3;
    % Evluation time step
    t = timeVector(nt) ;
    
    % For rigid systems these transformatio matrix system don`t change over time (this can be modify for system with yaw misalignment or tower deflection)
    r1   = expon( [thetaYaw ,          0, 0] );
    r2   = expon( [0        , -thetaTilt, 0] );
    r3   = eye(3,3);
    R12  = (r1*r2)*r3;  % Express local ref nacell coordinate into global coordinate
    R21  = R12' ;
    R34  = expon( [0,  thetaCone, 0]); % Express local blade coordinates into rotor coordinates
    R43  = R34' ;

    Hhub = Xtower + R12*Zhub ; % Hub height in global coordiantes

    %For each element update induced wind and calculate aerodynamic loads
    for i=1:3
        % Time dependecy tranformation matrix
        R23 = expon( [0  ,  0, thetaWing(nt, i) ] );
        R14 =  R12 * R23 * R34;  % Express a local coordinate of the blade into global coordinate.
        R41 = R14'            ;  % Convert a global coordinate into local blade coordinate
        
        for j=1:ni

            % Update global position of the blade node rb(i)
            ID              = sectionID(j)       ;                  % Blade section ID
            ch              = chord(j)           ;                  % Blade section chord
            beta            = twist(j)           ;                  % Blade section twist
            rblade          = [radius(j), 0, 0]' ;                  % Blade section position vector in local system
            rG(:, j, i, nt) = ( R14*rblade + R12*Zhub + Xtower )';      % Blade section position vector in global system
            rRot            = R34*rblade ;                         % Blade section position vector in system 3
    
            %Estimate relative wind speed (velocity triangle)            
            if Share
                V0g   = shareWindProfile(t, rG(1, j, i, nt), Vhub, Hhub(1), nPot); % wind velocity in global coordinates share wind
            elseif Uniform
                V0g   = uniformWindProfile( t, rG(1, j, i, nt), Vhub ) ;            % wind velocity in global coordinates uniform wind
            elseif Turbsim
                % Interpolate coordinates velocity for actual section location
                V0(1, j, i, nt) = interp2( Xgrid, Ygrid, Vx(:,:, nt), rG(1, j, i, nt) , rG(2, j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0(2, j, i, nt) = interp2( Xgrid, Ygrid, Vy(:,:, nt), rG(1, j, i, nt) , rG(2, j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0(3, j, i, nt) = interp2( Xgrid, Ygrid, Vz(:,:, nt), rG(1, j, i, nt) , rG(2, j, i, nt) , 'linear' ); % Puedes ajustar el método de interpolación
                V0g           = [ V0(1, j, i, nt); V0(2, j, i, nt); V0(3, j, i, nt) ] ;
                Egl           = [ 0 0 1 ; 0 -1 0 ; 1 0 0];
                V0g           = Egl*V0g ;  % convert from turbsim to BEM coordinates
                if j == 7
                    if i == 1
                        Vxblade(nt) = V0(1, j, i, nt);
                        Vyblade(nt) = V0(2, j, i, nt);
                        Vzblade(nt) = V0(3, j, i, nt);
                    end
                end
            end

            V0l      = R41*V0g               ;  % Convert wind velocity in global coordinate to local blade coordinates                                      
            L2       = [0 0 0; 0 1 0; 0 0 1] ;
            L3       = expon( [-pi/2 0 0] )  ;
            
%             if j == 10 && i == 2 || j == 10 && i == 3
%                 disp('Check Point')
%             end

            W     = Wake( 1:3, j, i,  nt-1)          ;
            vRot  = R14*[0  omega(nt,1)*rRot(1)  0]' ;   % V rot element expressed in global coordinates
            [vRel, vRelPerp, vRelG] = computeBEMVpiRels( V0g, vRot, W, R41, L2, L3 );
            vRelTest    = L2*(V0l + R41*W - vRot)          ;
            thetaTwist  = ( thetaPitch(nt, i) + beta ) ;
            chref       = expon( [-deg2rad(thetaTwist), 0, 0] )*[0, -ch, 0]'  ;

            [ AoArad, flowAngle ]   = computeLocalFlowAngle(chref, vRel, L3 ) ;
   
            vrelNorm = norm(vRel) ;  
            AoA(j, i, nt)           = rad2deg(AoArad) ; % Convert to degress to interpolate in airfoils coefs    
            
            %interpolate to find lift and drag (static)      
            %interpolate the values of lift and drag to the different thicknesses                
                       
            clstat    =  interp1( aoast(:,ID), liftCoef(:,ID),  AoA(j, i, nt), 'linear', 'extrap' );
            cdstat    =  interp1( aoast(:,ID), dragCoef(:,ID),  AoA(j, i, nt), 'linear', 'extrap' );
            
            % Apply Stall Delay Model 
            % Apply Dynamic Stall Model
            if DSM
                if OyeDSM
                    [ clift(j, i, nt), fs(j, i, nt) ] = dynamicStallModel(fst, fs(j, i, nt-1), clinvisc, ...
                                                                        clfs, delta_t, vrelNorm, ch);
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
        
                        [ BLcoefs(:, j, i, nt) ] = DSMbeddoessLeishman(  AoA(j, i, nt), vrelNorm, delta_t, ch, ...
                                                          BLcoefs(:, j, i, nt-1), AoA0,  Cd0, CnAo, ...
                                                          clstat, cdstat, Cn1, Cn2, AoA1, AoA2, T , ...
                                                          boolLeishmann );
                        
                        clift(j, i, nt)  = BLcoefs(1, j, i, nt) ;
                        cdrag(j, i, nt)  = cdstat;
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
            lift  =  1/2 * rho * clift(j, i, nt) * norm(chord(j)) * vrelNorm * vRelPerp ;
            drag  =  1/2 * rho * cdrag(j, i, nt) * norm(chord(j)) * vrelNorm * vRel     ;

            liftTest = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*clift(j, i, nt); % Lift force
            dragTest = 0.5*rho*( ( vrelNorm )^2 )*chord(j)*cdrag(j, i, nt); % Drag force
            
            BEMforces =  uBEMnormtanForces((lift + drag) , flowAngle ) ;
            
            fld(1:3, j, i, nt) = lift + drag;
            py = BEMforces(3) ;
            pz = BEMforces(2) ;

            % For time step update induced wind (relaxation method)
            % Apply glauert heavy induced factor correction
            n  = R43*[ 0, 0, -1 ]';
            Vprime = norm( vRel + n*dot( Wake(1:3, j, i, nt-1), n ));
            a     = -Wake(3, j, i, nt-1)/vrelNorm    ;
            atest = ( vrelNorm - Vprime ) / vrelNorm ;
            if(a<0.333)
                fglau=1;
            else
                fglau=0.25*(5-3*a);
            end
            
            %% Apply prandtl tip and hub correction factor
            F = prandtlFactor(radius(j), R, Rhub, flowAngle);

            %% Update induced wind
            aux1   = Vprime;
            aux2   = aux1*4*pi*rho*radius(j)*F;
            
            WakeQS(1:3, j, i, nt)  = -3*lift/aux2; 
            
            if DWM
                [Wake(1:3, j, i, nt), WakeInt(1:3, j, i, nt)] = dynamicInflowModel(Wake(1:3, j, i, nt-1), WakeQS(1:3, j, i, nt-1), ...
                                                                                   WakeInt(1:3, j, i, nt-1), WakeQS(1:3, j, i, nt), ...
                                                                                   vrelNorm, radius(j), R, a, R41, delta_t);

                WakeQS(1:3, j, i, nt)  = R14*WakeQS(1:3, j, i, nt)  ;
                WakeInt(1:3, j, i, nt) = R14*WakeInt(1:3, j, i, nt) ;
                Wake(1:3, j, i, nt)    = R14*WakeQS(1:3, j, i, nt)  ;


            else
                WakeQS(1:3, j, i, nt) = R14*WakeQS(1:3, j, i, nt) ;
                Wake(1:3, j, i, nt)   = R14*WakeQS(1:3, j, i, nt) ;
            end
      
        end   %end blade section 
    end   %end number of blades

    % Calculate thrust and power  

    for nblade=1:3 
            aeroTrq(nt,nblade)  = bladeMomentum( radius(:) , fld(:, :,nblade, nt), R34 );  
            blThrust(nt,nblade) = bladeThrust( radius(:), fld(:, :,nblade, nt), R34 ); 
    end

    aeroPower(nt, 1)   = ( aeroTrq(nt, 1)  + aeroTrq(nt, 2)  + aeroTrq(nt, 3 ) )*omega(nt,1);        % aerodynamic power
    thrust(nt, 1)      = ( blThrust(nt, 1) + blThrust(nt, 2) + blThrust(nt, 3) )/1000;               % total thrust kN
    rotTrq(nt, 1)      = ( aeroTrq(nt, 1)  + aeroTrq(nt, 2)  + aeroTrq(nt, 3 ) );                    % total Aerodynmic torque Nm
    
    %Turn on control after dynamic models are working
    if controlON
            [ genTrq(nt, 1), lastTimeGC(nt, 1), genOmega(nt, 1) ] = genControl( nGen, omega(nt, 1), thetaPitch(nt, 1), ...
                               genTrq(nt-1, 1), genOmega(nt-1, 1), lastTimeGC(nt-1, 1), timeVector(nt) ) ;
            
            genPower(nt, 1)    = genTrq(nt, 1)*nGen*omega(nt,1)*genEff;
   
            if nt <= length(timeVector) - 1 && timeVector(nt) > 10
                [ thetaPitch(nt+1, :), lastTimePC(nt+1, 1), speedError(nt+1,1), integError(nt+1, 1) ] = pitchControl( nGen, omega(nt, 1), thetaPitch(nt, :), ...
                                                                                            integError(nt, 1), lastTimePC(nt, 1), timeVector(nt) ) ;     
            end
    end
    omega(nt+1, 1)  =  omega(nt, 1) + ( ( rotTrq(nt, 1) - (genTrq(nt, 1)*nGen) )/I ) * delta_t ;
end  %end iteration(timesim)

%
studyCase = 'turbSimWindVel_15ms';
startPlot = 1; 
tStep     = delta_t ;
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ; 
folderPathFigs = './figs/uBEM/Outputs/noBLM_turbSim_12ms' ;
mkdir(folderPathFigs) 

fig1 = figure(1);
labelTitle = ' Rotor power ' ;
hold on; grid on
plot(timeVector(startPlot:end-1), genPower(startPlot:end-1,1)', 'r-', 'LineWidth',lw,'markersize',ms );
hold on
plot(timeVector(startPlot:end-1), aeroPower(startPlot:end-1,1)', 'b-', 'LineWidth',lw,'markersize',ms );
hold on
plot(timeVector(startPlot:end-1), ones(1,length(timeVector(startPlot:end-1)))*5e6, 'k-', 'LineWidth', lw,'markersize',ms )
legend('Gen Power', 'Aero Power', 'Rated Power','location','Best')
labx=xlabel('Time (s)'); laby=ylabel('Power (W)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig1 = strcat(folderPathFigs, [ '/rotorPower_', studyCase,'.jpg' ]) ;
print(fig1, namefig1,'-dpng') ;

fig2 = figure(2);
labelTitle =  ' Aerodynamic Power' ;
hold on; grid on
plot(timeVector(startPlot:end-1), genTrq(startPlot:end-1)*nGen', 'r-', 'LineWidth',lw,'markersize',ms )
hold on
plot(timeVector(startPlot:end-1), rotTrq(startPlot:end-1)', 'b-', 'LineWidth',lw,'markersize',ms )
legend('Gen Momentum', 'Aero Momentum', 'location','Best')
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

% if Turbsim
%     fig6 = figure(6);
%     hold on; grid on
%     labelTitle = 'Velocidad de viento en hub';
%     plot(timeVector(startPlot:end-1), Vxhub(startPlot:end-1, 1), 'r--', 'LineWidth', lw,'markersize',ms )
%     hold on
%     plot(timeVector(startPlot:end-1), Vyhub(startPlot:end-1, 1), 'b--', 'LineWidth', lw,'markersize',ms )
%     hold on
%     plot(timeVector(startPlot:end-1), Vzhub(startPlot:end-1, 1), 'k--', 'LineWidth', lw,'markersize',ms )
%     labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad (m/s) ');
%     set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
%     set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
%     title(labelTitle)
%     namefig6 = strcat(folderPathFigs, ['/hubVelocity=', studyCase,'.jpg'] ) ;
%     print(fig6, namefig6,'-dpng') ;
%     
%     fig7 = figure(7);
%     hold on; grid on
%     labelTitle = 'Velocidad de viento en punta de pala';
%     plot(timeVector(startPlot:end-1), Vxblade(startPlot:end-1, 1), 'r--', 'LineWidth', lw,'markersize',ms )
%     hold on
%     plot(timeVector(startPlot:end-1), Vyblade(startPlot:end-1, 1), 'b--', 'LineWidth', lw,'markersize',ms )
%     hold on
%     plot(timeVector(startPlot:end-1), Vzblade(startPlot:end-1, 1), 'k--', 'LineWidth', lw,'markersize',ms )
%     labx=xlabel( 'Time (s)' ); laby=ylabel(' Velocidad (m/s) ');
%     set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
%     set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
%     title(labelTitle)
%     namefig7 = strcat(folderPathFigs, ['/bladeVelocity=', studyCase,'.jpg'] ) ;
%     print(fig7, namefig7,'-dpng') ;
% end

fig8 = figure(8);
hold on; grid on
labelTitle = 'Fuerza de empuje sobre el rotor';
plot(timeVector(startPlot:end-1), thrust(startPlot:end-1, 1)', 'r-', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' Fuerza (N) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig8 = strcat(folderPathFigs, ['/totalThrust=', studyCase,'.jpg'] ) ;
print(fig8, namefig8,'-dpng') ;

cl     = zeros(length(timeVector(1:end-1)), 2) ;
cd     = zeros(length(timeVector(1:end-1)), 2) ;
AoAout = zeros(length(timeVector(1:end-1)), 2) ;
WakeG  = zeros(3, length(timeVector(1:end-1))) ;

for i = 1:length(timeVector) - 1
    WakeG(:,i)    = Wake(:, 15, 1, i) ;
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
labelTitle = 'Drag coef';
plot(timeVector(1000:end-1), cd(1000:end,1), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1000:end-1), cd(1000:end,2), 'b ', 'LineWidth', lw,'markersize',ms )
labx=xlabel( 'Time (s)' ); laby=ylabel(' C_d ');
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig10 = strcat(folderPathFigs, ['/dragCoef=', studyCase,'.jpg'] ) ;
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
plot(timeVector(1000:end-1), AoAout(1000:end, 1), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1000:end-1), AoAout(1000:end, 2), 'b', 'LineWidth', lw,'markersize',ms )
legend('Radius - 24.05 m', 'Radius - 56.167 m', 'location','Best')
labx=xlabel('Time (s)'); laby=ylabel('AoA (deg)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig12 = strcat(folderPathFigs, ['/AoA=', studyCase,'.jpg'] ) ;
print(fig12, namefig12,'-dpng') ;

fig13 = figure(13);
hold on; grid on
labelTitle = 'Wake component blade 1 section 15';
plot(timeVector(1:end-1), WakeG(1, :), 'r', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1:end-1), WakeG(2, :), 'b', 'LineWidth', lw,'markersize',ms )
hold on
plot(timeVector(1:end-1), WakeG(3, :), 'k', 'LineWidth', lw,'markersize',ms )
legend('Wx', 'Wy', 'Wz', 'location','Best')
labx=xlabel('Time (s)'); laby=ylabel('Wake (m/s)');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig13 = strcat(folderPathFigs, ['/Wake=', studyCase,'.jpg'] ) ;
print(fig13, namefig13,'-dpng') ;

clear fig1; clear fig2; clear fig3; clear fig4; clear fig5; clear fig6; clear fig7; clear fig8; clear fig9; clear fig10; clear fig11; clear fig12; clear fig13;
