function setUserParamLiveImagingAN
% these are good for the 60X imaging
% cell nuc area for the 40X imaging is ~ 1000-2000
global userParam;

userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
userParam.colonygrouping = 130; 

userParam.areanuclow = 1100;
userParam.areanuchi = 9000;
userParam.areanuclow2 = 1000;%800

userParam.flag = 1;% two parameters below are determined dynamically
%userParam.areanuclow_unmerge = 3000 ;  % min area of the merged object to start the split ( very specific to the image, need to generalize(select this parameter based on each inpit image)
userParam.minnucfragment =800;       % 1800 make this relatively small to cut off the junk stuff, then after unmerge use areanuclow2
userParam.linedil = 5;                 % size of the strel to dilate the line cut
userParam.tocut = 125;% 2450 if less then this parameter, don't cut


userParam.probthresh_nuc = 0.75;
userParam.probthresh_cyto = 0.9;

userParam.dilate_cyto = 5;
userParam.erode_nuc = 8;% 10

userParam.areacytolow = 50;