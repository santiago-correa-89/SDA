function R = expon( t );

  R = eye(3) ;

  tMod = norm( t ) ;

  if tMod > 0
    Rsk = skew(t) ;
    R = R + sin( tMod ) / tMod * Rsk  +  2 * ( sin( tMod/2 ) / tMod )^2 * Rsk^2 ;
  end

end