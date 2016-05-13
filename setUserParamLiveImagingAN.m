function setUserParamLiveImagingAN

global userParam;

userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
userParam.colonygrouping = 130;

userParam.areanuclow = 600;
userParam.areanuchi = 9000;
userParam.areanuclow2 = 600;%800

userParam.flag = 1;% two parameters below are determined dynamically
%userParam.areanuclow_unmerge = 6800 ;  % min area of the merged object to start the split ( very specific to the image, need to generalize(select this parameter based on each inpit image)
%userParam.minnucfragment = 1000;       % 2000
userParam.linedil = 7;                 % size of the strel to dilate the line cut

userParam.probthresh_nuc = 0.98;
userParam.probthresh_cyto = 0.85;

userParam.dilate_cyto = 5;
userParam.erode_nuc = 8;% 10

userParam.areacytolow = 200;