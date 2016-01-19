function newTrackParamAN

global userParam;


userParam.L = 24;
userParam.sizeImg = [1024, 1024];

userParam.verboseCellTrackerEDS = 0;
userParam.minTrajLen = 9;%4
userParam.mergeGap = 10;%2 6
userParam.sclDstCost = [1 2];
userParam.minlength = 20;
userParam.mincyto = 0;
userParam.splineparam = 0.9;