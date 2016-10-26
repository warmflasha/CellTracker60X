%% get max projetions from the live cell z-images
direc = ('');
direc2 = ('');

positions = (0:24);
tg = 0;
chan = 2;
Nchoose = [];
for k=2:size(positions,2)
pos = positions(k);
 MaxProjTimeGroupsAN(direc,direc2,pos,tg,chan,Nchoose);
 
end


%%
parpool(4);
%%
direc1 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/nuclear';
direc2 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/cyto';
% direc1 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-10-20-pluriset/nuc_projections/';
% direc2 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-10-20-pluriset/cyto_projections/';
% direc1 ='/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/nuc_projections';
% direc2 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/cyto_projections';

ff1 = readAndorDirectory(direc1);% nuc
ff2 = readAndorDirectory(direc2);% cyto
discardarea =800; % 1200 for the 60X data  900 for the 40X
mag = 40;
chan = [0,1];
cellIntensity = 300;% for pluri   280
cellIntensity1 = 200;% for cyto pluri 180 % 300 for tiling 1
parfor ii = 40%2:length(ff1.p) 
    disp(['Movie ' int2str(ff1.p(ii))]);
    nucmoviefile = getAndorFileName(ff1,ff1.p(ii),0,0,0);   % nuc channel ( last function argument)
    fmoviefile = getAndorFileName(ff2,ff2.p(ii),0,0,1);     % cyto channel
    [~, cmask, nuc_p, fimg_p] = simpleSegmentationLoop(nucmoviefile,fmoviefile,mag,cellIntensity,cellIntensity1);
    stmp = strsplit(nucmoviefile,direc1);
    ifile = ['/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/nuclear/' stmp{end}(2:(end-4)) '_{simplesegm}.h5'];%stmp{end-1}(2:end)
    nmask = readIlastikFile(ifile);
    nmask = cleanIlastikMasks(nmask,discardarea);% area filter last argument
    [newmasks, colonies] = statsArrayToSplitMasks(nmask,nuc_p,fimg_p,cmask);
    outfile = [int2str(ff1.p(ii)) '_tile1BGan.mat'];
    saveLiveCellData(outfile,newmasks,cmask,colonies);
end
%% plot cell trajectories
trajmin = 3;%10

outfile = '39_tile1BGan.mat';% good positions, Feb set[0,1,3,4,5,7,8,9,10,12,14,18,21,22,24,26,27,28,29,31,32,33];

plotcelltraces(outfile,trajmin)








 