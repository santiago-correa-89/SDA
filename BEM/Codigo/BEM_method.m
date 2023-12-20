function [ A , Aprime , PHI , ALPHA, CL, CD, ERR , j ] = BEM_method (B, c , r , t, lambda, theta_t, theta_p, Rrot, Rhub, Tol , CD_CL , flag_CLCD, factorCorreccion, aoa, clift, cdrag, cmomt)

%% Inicializo parámetros 
a0       = 1/3 ;
a_prime0 = 0   ;
eps      = 1e-3;
erraux   = Tol + eps ;
lambda_r = lambda*r/Rrot;

%     file = defineFile(t, path) ;
%     %% Guardar datos de perfil
%     datos = importAirfoilData(path, file);

%% Cálculo solidez
sigma    = B*c/( 2*pi*r ) ;   

%% Inicio iteración
j          = 1 ;
a(j)       = a0 ;
aprime(j)  = a_prime0 ;

while erraux > Tol
    
    %% Determino ángulod de velocidad relativa y ángulo de ataque
    phi(j)   = atan( (1 - a(j)) / ((( 1 + aprime(j)) * lambda_r ) + eps) ) ;%Angulo de velocidad relativa ( Radianes )
    thetaTwist = deg2rad(theta_t) + theta_p ;
    alpha(j) = phi(j) - thetaTwist          ;% Angulo de ataque 
    
    %% Defino coeficiente de sustentación y de arrastre según ángulo de ataque
    if ( flag_CLCD == 0 ) % en caso de contar con aproximación para párametros ingresarlos manualmente
        Cl(j) = 0.5 ; % Aproximado en grados
        Cd(j) = 0.01;
    elseif ( flag_CLCD == 1 )
        [Cl(j), Cd(j), Cm(j)] = airfoilCoef( alpha(j), [aoa(:,t), clift(:,t), cdrag(:,t), cmomt(:,t)] );
    else
        disp ('Wrong choices') ;
    end
    
    %% Determino factor de carga tangencial y normal adimiensional
    Cn         = Cl(j) * cos(phi(j)) + Cd(j) * sin(phi(j)) ;
    Ct         = Cl(j) * sin(phi(j)) - Cd(j) * cos(phi(j)) ;
    
    %% Factor de Correción de Prandtl
    fTip = (B/2)*( ( Rrot - r ) / ( r * sin(phi(j))) ) ;
    fHub = (B/2)*( ( r - Rhub ) / ( r * sin(phi(j))) ) ;
    FTip = (2/pi)*acos( exp( -fTip ) );
    FHub = (2/pi)*acos( exp( -fHub ) );
    F    = FTip*FHub;
    
    %% Determino a y aprime a partir de ecaución de método BEM aplicando factor de correción de Prandtl
    
    if factorCorreccion == 0
        % Proceso cuando no se aplica corrección por factor de inducción
        % axial elevado
        aa(j)      = 1/ ( (4*F*sin(phi(j))^2 / (sigma*Cn)) + 1) ; 
        aaprime(j) = 1/ ( (4*F*sin(phi(j))*cos(phi(j)) / (sigma*Ct)) - 1) ;

    elseif factorCorreccion == 1
        % Proceso para el factor de corrección de Glauert
        K  = 4*F*sin(phi(j))^2 / (sigma*Cn);
        aa(j) = 1/( K + 1) ;
        if aa(j) > 1/3
            aa(j) = glauertSolution(K, aa(j));
        end
        aaprime(j) = 1/ ( ((4*F*sin(phi(j))*cos(phi(j)))/(sigma*Ct)) - 1 ) ;

    elseif factorCorreccion == 2
        ac = 0.2; 
        K  = 4*F*sin(phi(j))^2 / (sigma*Cn);
        % Proceso para el factor de corrección de Wilson and Walker 1984
        aa(j) = 1/( K + 1) ;
        if aa(j) > ac
            aa(j) = 1 + (K/2)*(1 - 2*ac) - (1/2)*sqrt((K*(1-2*ac) + 2)^2 + 4*(K*(ac)^2 - 1)) ;
        end
        aaprime(j) = 1/ ((4*F*sin(phi(j))*cos(phi(j))/(sigma*Ct)) - 1) ;   
    else
        % Manejo de caso incorrecto o desconocido
        error('Valor de factorCorreccion no válido.');
    end

    %% Cálculo de error entre pasos de iteración
    err(j)     = max( [ (abs( aa(j) - a(j) )/((a(j)))) , (abs( aaprime(j) - aprime (j))/(aprime(j))) ] );
    
    %% Guardo valor de error en vector auxiliar
    erraux      = err(j) ;
    
    if erraux > Tol
        %% Incremento paso de iteración
        j = j + 1 ;
        %% Actualizo valor de a y aprime
        
        a(j)      = ( aa (j-1) +  a(j-1) ) / 2;
        aprime(j) = ( aaprime (j-1) + aprime(j-1) ) / 2;
        
    end

end

%% Variables de salida del método
A       = a(j)     ;
Aprime  = aprime(j);
PHI     = phi(j)   ;
ALPHA   = alpha(j) ; 
CL      = Cl(j)    ;
CD      = Cd(j)    ;
ERR     = err(j)   ;

end