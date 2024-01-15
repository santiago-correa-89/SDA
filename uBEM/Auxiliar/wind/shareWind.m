function Wind = shareWind(t, y, Vhub, Rhub, n)

  Wind = [Vhub*( (y/Rhub)^(n) ), 0, 0]' ;
  
end