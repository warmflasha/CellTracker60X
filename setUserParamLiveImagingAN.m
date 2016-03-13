function setUserParamLiveImagingAN

global userParam;

userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
userParam.colonygrouping = 130;

userParam.areanuclow = 600;
userParam.areanuchi = 8000;
userParam.areanuclow2 = 600;%800

userParam.flag = 1;% 
userParam.areanuclow_unmerge =  5800 ;%5600 6500 4820;
userParam.areanuchigh_unmerge =  6100 ;

userParam.probthresh_nuc = 0.98;
userParam.probthresh_cyto = 0.95;

userParam.dilate_cyto = 5;
userParam.erode_nuc = 8;% 10

userParam.areacytolow = 200;