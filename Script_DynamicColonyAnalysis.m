% add the loop over positions
% save group of every n files as a separate position (n - umber of time
% groups)

zplane = [];

chan = [0 1];


paramfile = 'setUserParamLiveImagingAN';

ilastikDirec1 = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_pluriW0');
ilastikDirec2 = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_pluriW1');
imgDirec1 = ('/Users/warmflashlab/Desktop/MaxProjectionsPluri_42hrNov29/Projections_NuclearChannel');% already max projections
imgDirec2 = ('/Users/warmflashlab/Desktop/MaxProjectionsPluri_42hrNov29/Projections_CytoChannel');% already max projections
fr_stim = [];% for the pluri dataset is empty
[nums, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec2,'CytoMask');%make two ilastik directories
timegroups = 3;% 4 for the diff dataset(nov12) and three for the pluri dataset
positions = length(nums)/timegroups; % number of separate position numbers (start from 0)
positions = 0:(positions-1);% vector with position numbers
%positions = 1; 
for kk=1: length(positions)
    
    pos = positions(kk);
    outfile = 'Pluri_42hrs.mat';% basic name for all positions
peaks = nucCytoIlastik2peaksLoop(ilastikDirec1,ilastikDirec2,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);
end

 addShiftToPeaks(outfile,fr_stim); 


% peaks{fr_stim} = peaks{fr_stim-1};%make peaks{38} = peaks{37} since frame 38 are completely different cells there (at least for the Nov12 imaging set)
% save(outfile,'peaks','imgfiles','imgfilescyto','-append');


runTracker(outfile,'newTrackParam');
global userParam;
userParam.colonygrouping = 120;
% look at colonies around the stimulation frame (window of couple hours)?
cellsToDynColonies(outfile);


fldat = [2 3];
delta_t = 5; % 5 in minutes (for the differentiated dataset), 11 mins for the pluripotent condition)
p = fr_stim*delta_t/60;
colSZ = 1;
flag = 1;
resptime = 24;% in frames ( converted to hours later)
 % add the loop oved positions and colony sizes here
GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);
datafin = GetDynamicColonyStats(outfile,fr_stim,delta_t,flag,colSZ,resptime);
