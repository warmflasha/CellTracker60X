%%
% add the loop over positions
% save group of every n files as a separate position (n - umber of time
% groups)
zplane = [];
chan = [0 1];
paramfile = 'setUserParamLiveImagingAN';

ilastikDirec1 = ('/Users/warmflashlab/Desktop/Jan8setIlastikMasks_headless_DiffW0');
ilastikDirec2 = ('/Users/warmflashlab/Desktop/Jan8setIlastikMasks_headless_DiffW1');
imgDirec1 =('/Users/warmflashlab/Desktop/MaxProjections_Dif_Jan8run/W0');% already max projections
imgDirec2 =('/Users/warmflashlab/Desktop/MaxProjections_Dif_Jan8run/W1') ;% already max projections


[nums, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec2,'CytoMask');%make two ilastik directories

timegroups = 3;% 

positions = length(nums)/timegroups; % number of separate position numbers (start from 0)

positions = 0:(positions-1);% vector with position numbers

for kk=1: length(positions)
    
    pos = positions(kk);
    outfile = 'jan8set_newparams.mat';% basic name for all positions
peaks = nucCytoIlastik2peaksLoop(ilastikDirec1,ilastikDirec2,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);
end
%%
 %addShiftToPeaks(outfile,fr_stim); 


% peaks{fr_stim} = peaks{fr_stim-1};%make peaks{38} = peaks{37} since frame 38 are completely different cells there (at least for the Nov12 imaging set)
% save(outfile,'peaks','imgfiles','imgfilescyto','-append');


addShiftToPeaks(outfile,fr_stim); %don't need to run the shift for the
peaks{fr_stim} = peaks{fr_stim-1};
save(outfile,'peaks','imgfiles','imgfilescyto');


