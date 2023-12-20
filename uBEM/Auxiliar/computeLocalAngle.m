function theta = computeLocalAngle( thElem, thNod )

if abs(thElem) < 1e-5
    thElem = 0;
end

if abs(thNod) < 1e-5
    thNod = 0;
end

if thNod == 0
    theta = 0;
elseif thNod ~= 0
    if sign(thElem) == sign(thNod)
        if abs(thElem) <= abs(thNod)
            theta = sign(thNod)*abs(thNod-thElem) ; 
        elseif abs(thNod) < abs(thElem) 
            theta = sign(thElem)*abs(thElem-thNod) ;    
        end
    elseif sign(thElem) ~= sign(thNod)
        if thNod < 0
            theta = sign(thNod) * ( abs(thNod) + abs(thElem));
        elseif thNod > 0
            theta = sign(thNod) * ( abs(thNod) + abs(thElem));
        elseif thNod == 0
            theta = thElem ;
        end
    end
end