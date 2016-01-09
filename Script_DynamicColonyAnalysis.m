% add the loop over positions
% save every 4th file as a separate position

zplane = [];

chan = [0 1];
[nums, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec1,'CytoMask');%make two ilastik directories
timegroups = 4;
positions = length(nums)/timegroups; % number of separate position numbers (start from 0)
positions = 0:positions;% vector with position numbers
paramfile = 'setUserParamLiveImagingAN';
outfile = '10ngml22hrs.mat';% basic name for all positions
ilastikDirec1 = ('/Users/warmflashlab/Desktop/IlastikMasks_headlessW0');
ilastikDirec2 = ('/Users/warmflashlab/Desktop/IlastikMasks_headlessW1');
imgDirec1 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W0');% already max projections
imgDirec2 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W1');% already max projections
for kk=1: length(positions)
    pos = positions(kk);
 
peaks = nucCytoIlastik2peaksLoop(ilastikDirec1,ilastikDirec2,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);
addShiftToPeaks(outfile,fr_stim);
runTracker(outfile,'newTrackParam');
global userParam;
userParam.colonygrouping = 120;
cellsToDynColonies(outfile);
end


fr_stim = 38;
fldat = [2 3];
delta_t = 5; 
p = fr_stim*delta_t/60;
colSZ = 1;
flag = 1;
 % add the loop oved positions here
GetDynamicColonyTraces(matfile,fr_stim,fldat,delta_t,colSZ);
datafin = GetDynamicColonyStats(matfile,fr_stim,delta_t,flag,colSZ);
