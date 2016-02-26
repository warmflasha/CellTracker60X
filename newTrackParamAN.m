function newTrackParamAN

global userParam;


userParam.L =100 ;% cell diameter is around 80 ,assume that cells move on average about ~ 1.2 their diameter 
userParam.sizeImg = [1024, 1024];

userParam.verboseCellTrackerEDS = 0;
userParam.minTrajLen = 4;%4
userParam.mergeGap = 4;%2 
userParam.sclDstCost = [1 2];
userParam.minlength = 20;
userParam.mincyto = 0;
userParam.splineparam = 0.9;