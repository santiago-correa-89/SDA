function [ AoA, flowAngle ] = computeLocalFlowAngle(chord, V, L3)

td = V/norm(V);
tl = L3*td    ;

tch = chord/norm(chord)   ;

% Define local system aligned with the reference system of the element
ex  =  [ 1, 0 ,0 ];
ey  =  [ 0, 1 ,0 ];
ez  =  [ 0, 0 ,1 ];

% Compute flow angle
cosFlow   = dot( td, tch ) / ( norm(td) * norm(tch) );
sinFlow   = dot(cross( td, tch ), ex ) / ( norm( td ) * norm( tch ));
AoA       = sign( sinFlow ) * acos( cosFlow );

cosFlowAngle   = dot( td, -ey ) / ( norm(ey) * norm(td) );
sinFlowAngle   = dot(cross( td, -ey ), ex ) / ( norm( td ) * norm( ey ));
flowAngle      = sign( sinFlowAngle ) * acos( cosFlowAngle );

end
