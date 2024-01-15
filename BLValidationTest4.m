%% Beddoess and Leishmann Validation test applying experimental data available through NREL from wind
%tunnel tests conducted at The Ohio State University (OSU) Aeronautical and
%Astronautical Research Laboratory (Ramsay et al 1995).

% OSU modelling parameters S809 Airfoil

%Reynolds 1 mill    -> f = 0.60 Hz Clean mean 14 - 10 amplitud sine wave (file C10|100_s809.txt)
%Reynolds 1 mill    -> f = 0.62 Hz LEGR mean 14 - 10 amplitud sine wave (file G10|100_s809.txt)

clear all; close all; clc;

%% Add path needed for code
addpath( genpath( [ pwd '/uBEM'] ) );
addpath( genpath( [ pwd '/TestData'] ) );

%% S809 BL model inputs
boolLeishmann = true ;
Chord         = 0.457    ;  % Airfoil chord (m)
Cd0           = 0.0017   ;  % Zero lift drag coef
AoA0          = 0.908    ;  % zero lift angle of attack
AoA1          = 8.57     ;  % critical stall angle AoA > AoA0  
AoA2          = 18.1     ;  % critical stall angle AoA < AoA0
CnSlope       = 6.684    ;  % Cn AoA slope of linear Cl region (1/rad)
Cn1           = 1.0093   ;  % Critical value of Cn′ at LE sepearation , Cn1. It can be calculated   
                                  % static value of Cn at either the break in the pitching
                                  % from the moment or the loss of chord force at the onset of stall
Cn2           = 2.4500   ;  % Critical value of Cn′ at LE separation for negative AOAs

%% Initialization parameters
state     = 'clean' ;
testAmp   = '10'; 
testMean  = '14'; 
% S809 airfoil data
staticS809dataRe1m        = table2array( readtable('/TestData/static_s809_Re1m.txt'                 ) );
OSUdata   = table2array( readtable('/TestData/unesteadyClean_Mean14_Amp10_Re1m.txt'  ) );
%OSUdata    = table2array( readtable('/TestData/unesteadyLEGR_Mean14_Amp10_Re1m.txt'   ) );

%% Experimental conditions
vrelNorm = 109.9*0.3048           ;  % Wind speed (m/s) - 109.9 - 110.4
fHz      = 0.60                   ;  % Oscillator frequency (hz)
w        = fHz*2*pi               ;  % Oscillator frequency (rad/s)
k        = w*Chord/( 2*vrelNorm ) ;  % Reduce frequency

%% Static Experimental data Test OSU
AoAst = staticS809dataRe1m(:,1); 
Clst  = staticS809dataRe1m(:,2);
Cdst  = staticS809dataRe1m(:,3);
for i = 1:length(AoAst)
    Cnst(i)  = Clst(i)*cos(AoAst(i)*pi/180)  +  (Cdst(i) - Cd0)*sin(AoAst(i)*pi/180) ;
end

%% Unsteady Experimental data Test OSU
expT    = OSUdata(:,2) ;
expAoA  = OSUdata(:,3) ;
expCl   = OSUdata(:,4) ;
expCd   = OSUdata(:,5) ;

%% variables
tVector = round(0:( max(expT)/119 ):max(expT), 3 ); % Time vector
tIter   = length(tVector);
AoA     = 14 - 10*cos(w*tVector) ; % Angle of attack variation

%% Initiliazation BL model
nBlades     = 1;
nSections   = 1;
T           = [ 2.8, 3.5, 7.0, 9.0 ] ;    %  LBM model time constant Tp, Tf0, Tv0, Tvl ( ref Pereira 2011 )
[ BLcoefs ] = initDSMbeddoesLeishman( nSections, nBlades, tIter ) ;

for i = 2:tIter
    for j = 1:nBlades
        for k = 1:nSections
            AoAn   = AoA(i);
            deltaT = tVector(i) - tVector(i-1) ;
            clstat = interp1( AoAst(:), Clst(:),  AoAn, 'linear', 'extrap' );
            cdstat = interp1( AoAst(:), Cdst(:),  AoAn, 'linear', 'extrap' );

            [ BLcoefs(:, k, j, i) ] = DSMbeddoessLeishman(  AoA(i), vrelNorm, deltaT, Chord, ...
                                                      BLcoefs(:, k, j, i-1), AoA0,  Cd0, CnSlope, ...
                                                      clstat, cdstat, Cn1, Cn2, AoA1, AoA2, T , ...
                                                      boolLeishmann );
        end
    end
    clift(:, i)  = BLcoefs(1, :, :, i) ;
    cdrag(:, i)  = BLcoefs(2, :, :, i) ;
    cnorm(:, i)  = BLcoefs(3, :, :, i) ;
    cchord(:, i) = BLcoefs(4, :, :, i) ;
    CnP(:, i)    = BLcoefs(15,:, :, i) ;
    CnF(:, i)    = BLcoefs(22,:, :, i) ;
    CnV(:, i)    = BLcoefs(19,:, :, i) ;
    Cvn(:, i)    = BLcoefs(21,:, :, i) ;
end

errCl     = abs((expCl(60:end) - clift(60:end)')./expCl(60:end) )*100 ;
maxErrCl  = max(errCl')                ;
meanErrCl = mean(errCl')./expCl(60:end)        ;
mseErrCl  = immse(expCl(60:end), clift(60:end)')      ;
rmseErrCl = rmse(expCl(60:end), clift(60:end)')       ;

errCd     = abs(( expCd(60:end) - cdrag(60:end)')./expCd(60:end))*100 ;
maxErrCd  = max(errCd)                ;
meanErrCd = mean(errCd)./expCd(60:end)        ;
mseErrCd  = immse(expCd(60:end), cdrag(60:end)')       ;
rmseErrCd = rmse(expCd(60:end), cdrag(60:end)')        ;

lw = 1.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 1 ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ; 
folderPathFigs = ['./figs/LBmodel/validationTest/', state, 'amp', testAmp, 'mean', testMean] ;
mkdir(folderPathFigs) 

fig1 = figure(1);
hold on; grid on
labelTitle = 'Uneastedy C_l case with clean wind tunnel';
plot(AoAst(11:31), Clst(11:31), 'r', 'LineWidth', 2, 'markersize', ms )
hold on
plot(AoAst(11:31), CnSlope*pi/180*(AoAst(11:31) + AoA0), 'g--', 'LineWidth', lw, 'markersize', ms )
hold on
plot(expAoA(60:100), expCl(60:100), '-.k', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(50:end), clift(50:end), '-.m', 'LineWidth', lw,'markersize',ms )
legend('Steady Cl data', 'Steady Cn data', 'Cn slope','OSU data', 'B-L model','Best')
labx=xlabel( 'AoA (º)' ); laby=ylabel(' C_l ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig1 = strcat(folderPathFigs, ['/liftCoef_', state, state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig1, namefig1,'-dpng') ;

fig2 = figure(2);
hold on; grid on
labelTitle = 'Uneastedy C_d case with clean wind tunnel';
plot(AoAst(15:31), Cdst(15:31), 'r', 'LineWidth', lw, 'markersize', ms )
hold on
plot(expAoA(60:end), expCd(60:end), '--ob', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(60:end), cdrag(60:end), '--+k', 'LineWidth', lw,'markersize',ms )
legend('Steady Cd data', 'OSU data', 'B-L model','Best')
labx=xlabel( 'AoA (º)' ); laby=ylabel(' C_d ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig2 = strcat(folderPathFigs, ['/dragCoef_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig2, namefig2,'-dpng') ;

fig3 = figure(3);
hold on; grid on
labelTitle = 'Uneastedy C_l case with clean wind tunnel';
plot(expAoA(60:end), expCl(60:end), 'k', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(60:end), clift(60:end), '--+b', 'LineWidth', lw,'markersize',ms )
legend( 'OSU data', 'B-L model', 'Best')
labx=xlabel( 'AoA (º)' ); laby=ylabel(' C_l ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig3 = strcat(folderPathFigs, ['/versusliftCoef_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig3, namefig3,'-dpng') ;

fig4 = figure(4);
hold on; grid on
labelTitle = 'Uneastedy C_d case with clean wind tunnel';
plot(expAoA(60:end), expCd(60:end), 'k', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(60:end), cdrag(60:end), '-+b', 'LineWidth', lw,'markersize',ms )
legend('OSU data', 'B-L model', 'Best')
labx=xlabel( 'AoA (º)' ); laby=ylabel(' C_d ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig4 = strcat(folderPathFigs, ['/versusdragCoef_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig4, namefig4,'-dpng') ;

fig5 = figure(5);
hold on; grid on
labelTitle = 'AoA variation Experimental vs Model';
plot(tVector, AoA, '--or', 'LineWidth', lw, 'markersize', ms )
hold on
plot(expT, expAoA, '-+b', 'LineWidth', lw,'markersize',ms )
legend('AoA model', 'OSU exp','Best')
labx=xlabel( 'Time (s)' ); laby=ylabel(' AoA (º) ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig5 = strcat(folderPathFigs, ['/AoAvariation_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg']  ) ;
print(fig5, namefig5,'-dpng') ;

fig6 = figure(6);
hold on; grid on
labelTitle = 'Uneastedy C_l models result';
plot(AoAst(11:31), Clst(11:31), 'r', 'LineWidth', lw, 'markersize', ms )
hold on
plot(AoA(60:100), clift(60:100), 'k', 'LineWidth', lw, 'markersize', ms )
hold on
plot(AoA(60:100), CnP(60:100), 'b--', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(60:100), CnF(60:100), 'g--', 'LineWidth', lw,'markersize',ms )
hold on
plot(AoA(60:100), CnV(60:100), 'm--', 'LineWidth', lw,'markersize',ms )
legend('Steady Cl data', 'Total model', 'Attached Component', 'Viscuous component', 'Vortex advection component','Best')
labx=xlabel( 'AoA (º)' ); laby=ylabel(' C_l ');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig6 = strcat(folderPathFigs, ['/components_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig6, namefig6,'-dpng') ;

fig7 = figure(7);
hold on; grid on
labelTitle = 'C_l error';
edges = linspace(min(errCl), max(errCl), 20);
histogram(errCl, edges, 'FaceColor', [0.6 0.8 1]);
labx=xlabel( 'Error (%)' ); laby=ylabel(' Instances');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
title(labelTitle)
namefig7 = strcat(folderPathFigs, ['/errorCl_', state, '_Amp', testAmp, '_mean', testMean, '_Re1m.jpg'] ) ;
print(fig7, namefig7,'-dpng') ;