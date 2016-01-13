function [ncells, cellmrna, mrnapercell] = runmaskone(segfiledir, rawfiledir, nzslices, colonyno, objno)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices


[dirinfo, start] = readdirectory(segfiledir);
[dirinfo1, start1] = readdirectory(rawfiledir);


maskno = start + (colonyno-1)*nzslices; % ilastik mask
imageno = start1 + (colonyno-1)*nzslices ; %nuclear channel raw image

[pnuc, inuc] = readmaskfiles(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
%%
logfilter = 10;
bthreshfilter = 0.25;
diskfilter = 3;
area1filter = 100;
pmasks = primaryfilter(pnuc,logfilter, bthreshfilter, diskfilter, area1filter);
%%
minstartobj = 4;
minsolidity = [0.9, 0.8];
area2filter = 300;
[zrange, smasks] = secondaryfilter(pmasks, minstartobj, minsolidity, diskfilter, area2filter);
%%
zmatch = 12;
matchdistance = 15;
[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, matchdistance, zrange, zmatch);
%%
overlapthresh = 80;
colonyno = (maskno - start)/nzslices + 1;
[ncells, nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, overlapthresh);

%%
mrnamatfile = '/Volumes/data/Sapna/150813FISH_MP/200um/ch3allspots.mat';
[cellmrna, mrnapercell] = mrnapercells(nucleilist, stats, mrnamatfile, colonyno, zrange);

%%
if (~exist('objno'))
    objno = [5];
end
nucleimrnacheck(masterCC, inuc, nucleilist, stats, objno, cellmrna, zrange);
