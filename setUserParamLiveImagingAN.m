function setUserParamLiveImagingAN

global userParam;

userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
userParam.colonygrouping = 130;

userParam.areanuclow = 1200;
userParam.areanuchi = 15000;
userParam.areanuclow2 = 1200;

userParam.flag = 1;
userParam.areanuclow_unmerge = 6500;%4820;

userParam.probthresh_nuc = 0.98;
userParam.probthresh_cyto = 0.98;

userParam.dilate_cyto = 5;
userParam.erode_nuc = 10;

userParam.areacytolow = 50;