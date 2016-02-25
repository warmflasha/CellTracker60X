function runmaskonecheck(segfiledir, rawfiledir, paramfile, samplenum, position_num)
%%
tic;
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices
% samplenum: no. of different samples/conditions compared.
% position_num: within a particular sample, th position for which you want
% to check the segmentation.
% objno: cell no. for which you want to check no. of mRNA's assigned (for

global userparam
eval(paramfile);

ff = readFISHdir(rawfiledir, userparam.nsamples);
nzslices = ff.zslices{samplenum}(position_num);

nuc_ch = userparam.nucchannel;
[pnuc, inuc] = readmaskfilesnew(segfiledir, rawfiledir, samplenum, position_num,  nzslices, nuc_ch);
%%
pmasks = primaryfilter(pnuc, userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
%%

[zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
%%

[PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userparam.matchdistance, zrange, userparam.zmatch);
%%

[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);

toc;