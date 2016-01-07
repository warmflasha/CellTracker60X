% add the loop over positions
fr_stim = 38;
zplane = [];
pos = 1;
chan = [0 1];
fldat = [2 3];
delta_t = 5;
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
p = fr_stim*delta_t/60;
colors = colorcube(50);
for i = 1:size(colonies,2);
ratio = cell(size(colonies(i).cells,2),1);
tpt = cell(size(colonies(i).cells,2),1);

    for k=1:size(colonies(i).cells,2)
        ratio{k} = colonies(i).cells(k).fluorData(:,fldat(1))./colonies(i).cells(k).fluorData(:,fldat(2));
        tpt{k} =  (colonies(i).cells(k).onframes');
        tpt{k} =  (tpt{k}.*delta_t)./60;
    end
    for k=1:size(ratio,1)
       figure(i), plot(tpt{k},ratio{k},'-*','color',colors(k,:));hold on
        legend(['bmp4 added at ' num2str(p) 'hours']);
    end
    ylabel('nuc/cyto raio');
    xlabel('Time, hours');
end



