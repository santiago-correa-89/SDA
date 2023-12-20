function [angle , eVel ]= computeFlowAngle(V)

L3 = expon( [0 -pi/2 0] ) ;
ed = V/norm(V);
el = L3*ed    ;

% Define local system aligned with the reference system of the element
ex  =  [ 1, 0 ,0 ];
ey  =  [ 0, 1 ,0 ];
ez  =  [ 0, 0 ,1 ];

% Compute flow angle
cosFlow   = dot( ed, -ez ) / ( norm(ez) * norm(ed) );
sinFlow   = dot(cross( ed, -ez ), ey ) / ( norm( ed ) * norm( ez ));
angle     = sign( sinFlow ) * acos( cosFlow );

eVel = [el, ey', ed];

end
