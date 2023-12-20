function [ genTrq, lastTimeGC, genOmegaF ] = genControl( nGen, rotOmega, bladePitch, lastGenTrq, genOmegafTP, lastTimeGC, t1p )

deltaTime  = t1p - lastTimeGC;
genOmega   = nGen*rotOmega;  % actual generator speed


if deltaTime < 0.05

    genOmegaF  = genOmegafTP;
    lastTimeGC = lastTimeGC ;
    genTrq     = lastGenTrq ;

elseif deltaTime >= 0.05

    %% Filter the generator speed
    cornerFreq = 1.570796                                     ; % Frequency parameter of filter
    alpha      = exp( ( deltaTime )*cornerFreq )              ; % Constant parameter of filter
    genOmegaF  = genOmega ; % filter gen spee, genOmegafTP previous stepp

    %% Parameters
    omegaRatedSp  = 121.6805   ;    % Rated generator speed in rad/seg
    cutInSpeed    = 70.16224   ;    % Cut-in omega in rad/seg
    Rg15_2LimitSp = 91.21091   ;    % Region 1 1/2 upper limit speed rad/seg
    minPitch      = 0.01745329 ;    % 1 deg pitch angle compute the region 3 in that case
    ratedSynSp    = 10.0       ;    % Rated generator slip percentage in Region 2
    maxTorque     = 47402.91   ;    % Max generator torque in region 3
    ratedPower    = 5296610.0  ;    % Rated generator power in region 3
    maxRatedTorq  = 15000.0    ;    % Maximum torque rate
    rgn2const     = 2.332287   ;    % Generator torque constant in Region 2

    %% Torque parameters
    syncSpeed = omegaRatedSp/( 1.0 + 0.01*ratedSynSp ) ;
    slope15   = ( rgn2const * Rg15_2LimitSp * Rg15_2LimitSp )/(Rg15_2LimitSp - cutInSpeed ) ;
    slope25   = ( ratedPower/omegaRatedSp )/( omegaRatedSp - syncSpeed ) ;

    if ( rgn2const == 0.0 ) % if the Region 2 torque is flat, and thus, the denominator in the ELSE condition is
        Rg2LimitSp = syncSpeed ;
    else % if the Region 2 torque is quadratic with speed
        Rg2LimitSp = ( slope25 - sqrt( slope25*( slope25 - 4.0*rgn2const*syncSpeed ) ) )/( 2.0*rgn2const ) ;
    end

    %% Compute the generator torque, which depends on which region we are in:
        
        if genOmegaF >= omegaRatedSp || (bladePitch >= minPitch)% Region 3 - power is constant
            
            genTrq = ratedPower / genOmegaF;
        
        elseif genOmegaF <= cutInSpeed  % Region 1 - torque is zero
            
            genTrq = 0.0;
        
        elseif genOmegaF > cutInSpeed && genOmegaF <= Rg15_2LimitSp  % Region 1 1/2 - linear ramp in to
            
            genTrq = slope15 * (genOmegaF - cutInSpeed);
        
        elseif genOmegaF > Rg15_2LimitSp && genOmegaF <= Rg2LimitSp  % Region 2 - optimal torque is pro
            
            genTrq = rgn2const * genOmegaF * genOmegaF;
        
        elseif (genOmegaF > Rg2LimitSp) && (genOmegaF < omegaRatedSp)% We are in region 2 1/2 - simple induction
            
            genTrq = slope25 * (genOmegaF - syncSpeed);
        end

    
    % Saturate the commanded torque using the maximum torque limit:
    genTrq = min( genTrq, maxTorque );
    
    trqRate = ( genTrq - lastGenTrq ) / deltaTime;               % Torque rate (unsaturated)
    trqRate = min ( max( trqRate, -maxRatedTorq), maxRatedTorq); % Saturate the torque rate using its maximum absolute value
    genTrq  = lastGenTrq + trqRate * deltaTime;                  % Saturate the command using the torque rate limit
    % Reset the values of LastTimeVS and LastGenTrq to the current values:
    lastTimeGC      = t1p;
end
end
