function [aoa,cl,cd,cm,fstat,clinv,clfullysep]=...
    readairfoildata_v3(filename1,filename2,...
    filename3,filename4,filename5,filename6)
strucdata=importdata(filename1);
aoa(:,1)=strucdata(:,1);
cl(:,1)=strucdata(:,2);
cd(:,1)=strucdata(:,3);
cm(:,1)=strucdata(:,4);
fstat(:,1)=strucdata(:,5);
clinv(:,1)=strucdata(:,6);
clfullysep(:,1)=strucdata(:,7);

strucdata=importdata(filename2);
aoa(:,2)=strucdata(:,1);
cl(:,2)=strucdata(:,2);
cd(:,2)=strucdata(:,3);
cm(:,2)=strucdata(:,4);
fstat(:,2)=strucdata(:,5);
clinv(:,2)=strucdata(:,6);
clfullysep(:,2)=strucdata(:,7);

strucdata=importdata(filename3);
aoa(:,3)=strucdata(:,1);
cl(:,3)=strucdata(:,2);
cd(:,3)=strucdata(:,3);
cm(:,3)=strucdata(:,4);
fstat(:,3)=strucdata(:,5);
clinv(:,3)=strucdata(:,6);
clfullysep(:,3)=strucdata(:,7);

strucdata=importdata(filename4);
aoa(:,4)=strucdata(:,1);
cl(:,4)=strucdata(:,2);
cd(:,4)=strucdata(:,3);
cm(:,4)=strucdata(:,4);
fstat(:,4)=strucdata(:,5);
clinv(:,4)=strucdata(:,6);
clfullysep(:,4)=strucdata(:,7);

strucdata=importdata(filename5);
aoa(:,5)=strucdata(:,1);
cl(:,5)=strucdata(:,2);
cd(:,5)=strucdata(:,3);
cm(:,5)=strucdata(:,4);
fstat(:,5)=strucdata(:,5);
clinv(:,5)=strucdata(:,6);
clfullysep(:,5)=strucdata(:,7);

strucdata=importdata(filename6);
aoa(:,6)=strucdata(:,1);
cl(:,6)=strucdata(:,2);
cd(:,6)=strucdata(:,3);
cm(:,6)=strucdata(:,4);
fstat(:,6)=strucdata(:,5);
clinv(:,6)=strucdata(:,6);
clfullysep(:,6)=strucdata(:,7);