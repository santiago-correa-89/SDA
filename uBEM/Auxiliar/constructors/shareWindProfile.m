function Wind = shareWindProfile(t, x, Vhub, Rhub, n)

  Wind = [0, 0, Vhub*( (x/Rhub)^(n) )]' ;
  
end