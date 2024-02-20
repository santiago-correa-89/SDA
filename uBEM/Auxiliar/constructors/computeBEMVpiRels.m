function [VpiRel, VpiRelPerp, VrelG] = computeBEMVpiRels( udotFlow, udotFrame, udotWake, R0, L2, L3 )
  
% the relative velocity in global cooridantes is Ref: eq 9.9 Hansen 2015
VrelG = udotFlow - udotFrame + udotWake;

% then the projection (in t2,t3 plane) of the relative flow velocity in deformed coordinates is:
VpiRel = L2*R0*VrelG ;
  
% the perpendicular flow relative velocity projection in deformed coordinates is:
VpiRelPerp = L3 * VpiRel ;

end