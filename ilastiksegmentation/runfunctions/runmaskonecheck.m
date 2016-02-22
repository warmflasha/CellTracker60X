function runmaskonecheck(segfiledir, rawfiledir,  position_num, paramfile)
%%
tic;
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices
% position_num: image no (position_num=1, implies I image or 1 unique position imaged)
% objno: cell no. for which you want to check no. of mRNA's assigned (for

global userparam
eval(paramfile);

[dirinfo, start] = readdirectory(segfiledir);
[dirinfo1, start1] = readdirectory(rawfiledir);


maskno = start + (position_num-1)*userparam.nzslices; % ilastik mask
imageno = start1 + (position_num-1)*userparam.nzslices ; %nuclear channel raw image
[pnuc, inuc] = readmaskfiles(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, userparam.nzslices, imageno);

%%
pmasks = primaryfilter(pnuc, userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
%%

[zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
%%

[PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userparam.matchdistance, zrange, userparam.zmatch);
%%

[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imview);

toc;