function setUserParam3DsegmentationAN

global userParam;

% parameters for the 3d segmentation ( sapna/idse code)

userParam.logfilter = 10;
userParam.bthreshfilter = 0.25;
userParam.diskfilter = 4;
userParam.area1filter = 100;

userParam.minstartobj = 1;
userParam.minsolidity = [0.9, 0.8];
userParam.area2filter = 700;

userParam.zmatch = 5;
userParam.matchdistance = 25;

userParam.overlapthresh = 60;