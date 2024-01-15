function [uturb]=turbintpol(thetapos,NA,radpos,NR,thetapol,radpol,...
    ntime,upol);
    %Determine the azimuthal position of the blade, relative to polar grid
    azim=thetapos;
    azim=mod(azim,2*pi);
    for j=2:NA+1;
        if( (azim<=thetapol(j)) & (azim>=thetapol(j-1)))
           jpos=j;               
        end
    end
    
        
    %Azimuthal interpolation     
    for ri=1:NR+1;
      u1(ri,ntime)=upol(ri,jpos,ntime);
      u2(ri,ntime)=upol(ri,jpos-1,ntime);
    
      w1=thetapol(jpos)-azim;
      w2=azim-thetapol(jpos-1);
      norm=w1+w2;
      w1=w1/norm;
      w2=w2/norm;
      urad(ri,ntime)=w2*u1(ri,ntime)+w1*u2(ri,ntime);
    end
      
      %Radial interpolation
      uturb=interp1(radpol(:),urad(:,ntime),radpos);
