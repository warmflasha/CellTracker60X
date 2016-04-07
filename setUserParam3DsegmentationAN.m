function setUserParam3DsegmentationAN

global userParam;

% parameters for the 3d segmentation ( sapna/idse code)

userParam.logfilter = 9;
userParam.bthreshfilter = 0.25;% 0.25
userParam.diskfilter = 3;%3  4
userParam.area1filter = 600;

userParam.minstartobj = 1;
userParam.minsolidity = [0.9, 0.8];%[0.9 0.8]
userParam.area2filter = 800;%1000

%userParam.zmatch = 4;% thisparameter is set in the function to be the size
%of the zrange (so that the foud nuclei would be traced throughout all
%nonempty zplanes)
userParam.matchdistance = 20;%15 25

userParam.overlapthresh = 80;