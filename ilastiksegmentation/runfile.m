function runfile(segfiledir, rawfiledir, nzslices, outfile)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices


[dirinfo, start] = readdirectory(segfiledir);
[dirinfo1, start1] = readdirectory(rawfiledir);

for maskno = start:nzslices:size(dirinfo,1)
    colonyno = (maskno - start)/nzslices + 1;
    imageno = start1 + (colonyno-1)*nzslices ;
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
    [ncells(colonyno), nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, overlapthresh);
    
    
    %%
    
    for ch = 1:2
        mrnamatfile = sprintf('/Volumes/data/Sapna/150813FISH_MP/200um/ch%dallspots.mat', ch);
        [cellmrna{colonyno}{ch}, mrnapercell{colonyno}{ch}] = mrnapercells(nucleilist, stats, mrnamatfile, colonyno, zrange);
        
    end
end

save(outfile, 'cellmrna', 'mrnapercell', 'ncells');
end
