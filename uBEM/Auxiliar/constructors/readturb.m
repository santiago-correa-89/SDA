function [NTime,NA,NR,radpol,thetapol,delta_t,maxtime,upol]=...
    readturb(file)
fid=fopen(file,'r');
pnum=fread(fid,1,'double');
NTime=fread(fid,1,'double');
NA=fread(fid,1,'double');
NR=fread(fid,1,'double');
time=fread(fid,[NTime 1],'double');
radpol=fread(fid,[NR+1 1],'double');
thetapol=fread(fid,[NA 1],'double');
%time=fread(fid,[NTime 1],'double');
usim=fread(fid,[pnum,NTime],'double');
fclose(fid);



% fid=fopen('upol.bin','w');
% fwrite(fid,pnum,'double');
% fwrite(fid,NTime,'double');
% fwrite(fid,NA,'double');
% fwrite(fid,NR,'double');
% fwrite(fid,time,'double');
% fwrite(fid,radpol,'double');
% fwrite(fid,thetapol,'double');
% fwrite(fid,usim,'double');
% fclose(fid);










delta_t=time(2)-time(1)
maxtime=NTime

%Velocity put in a polar grid (radpol,thetapol)
for ntime=1:NTime
  for j=1:NA
    % First point is the center
    upol(1,j,ntime)=usim(1,ntime);
    % Axial velocity in polar grid (i,j)
    for i=2:NR+1
      upol(i,j,ntime)=usim(i+NR*(j-1),ntime);
    end
  end
end

%Ensuring that the last azimuthal element is the same as the first
for ntime=1:NTime;
    for i=1:NR+1
        upol(i,NA+1,ntime)=upol(i,1,ntime);
    end
end
thetapol(NA+1)=2*pi;
