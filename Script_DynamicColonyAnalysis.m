% add the loop over positions

zplane = [];
pos = 7;
chan = [0 1];

paramfile = 'setUserParamLiveImagingAN';
outfile = 'testheadless.mat';% basic name for all positions
ilastikDirec = ('/Users/warmflashlab/Desktop/TestTraces');
imgDirec1 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W0');% already max projections
imgDirec2 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W1');% already max projections

peaks = nucCytoIlastik2peaksLoop(ilastikDirec,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);
addShiftToPeaks(outfile,fr_stim);
runTracker(outfile,'newTrackParam');
global userParam;
userParam.colonygrouping = 120;
cellsToDynColonies(outfile);


%colonies(k).cells
%colonies(k).cells(j).fluorData
%plot(colonies(1).cells(2).fluorData(:,2)./colonies(1).cells(2).fluorData(:,3),'r*--');
%vect{j} = 1:length(peaks);
%vect{j} = (vect{j}.*delta_t)./60;
fr_stim = 38;
fldat = [2 3];
delta_t = 5; 
p = fr_stim*delta_t/60;
colSZ = 1;
flag = 1;
 % add the loop oved positions here
GetDynamicColonyTraces(matfile,fr_stim,fldat,delta_t,colSZ);
datafin = GetDynamicColonyStats(matfile,fr_stim,delta_t,flag,colSZ);
