function [ BLcoefs ] = initDSMbeddoesLeishman(nSections, nBlades, tIter)

    BLcoefs               = zeros(22, nSections, nBlades, tIter)  ;
    BLcoefs(16, :, :,  1) = 1                 ; % fnPrime intialized as 1

end