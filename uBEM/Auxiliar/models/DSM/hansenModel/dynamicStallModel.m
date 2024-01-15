function [Cl, fs] = dynamicStallModel(fst, fsold, clinvisc, clfs, deltat, Vrel, c)

A     = 4;
Tau0  = A*c/sqrt(Vrel);

fs    = fst + (fsold - fst) * exp(-deltat/Tau0);
        
% Dynamic stall lift correction
Cl = clinvisc * fs + clfs * (1 - fs) ;

end