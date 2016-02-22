function runmaskdir(segfiledir, rawfiledir, mrnafilepath, paramfile, outfile, fluorpdir)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices: no. of z slices
% outfile: final output file path

global userparam
eval(paramfile)

[dirinfo, start] = readdirectory(segfiledir);
[dirinfo1, start1] = readdirectory(rawfiledir);


nzslices = userparam.nzslices;

for maskno = start:nzslices:size(dirinfo,1)
    colonyno = (maskno - start)/nzslices + 1;
    imageno = start1 + (colonyno-1)*nzslices ;
    [pnuc, inuc] = readmaskfiles(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
    %%
    
    pmasks = primaryfilter(pnuc,userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
    %%
   
    [zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
    %%
    
    [PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userparam.matchdistance, zrange, userparam.zmatch);
    %%
    
    colonyno = (maskno - start)/nzslices + 1;
    
    [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);
        
    %%
    
    peaks1{colonyno} = mrnapercells(nucleilist, stats, mrnafilepath, colonyno, zrange, userparam.channels, userparam.cmcenter);
    %%
    peaks{colonyno} = peakscelltrackerformat(peaks1{colonyno});
    
    %%
    if(exist(fluorpdir) && colonyno ~= userparam.negativecontrol)
        [dirinfo2, start2] = readdirectory(fluorpdir);
        
        if colonyno < userparam.negativecontrol
            fimageno = start2+(colonyno-1)*nzslices;
        else
            fimageno = start2+(colonyno-2)*nzslices;
        end
        peaks{colonyno} = masksin3d(CC, nucleilist, nzslices, fimageno, fluorpdir, dirinfo2, inuc, peaks{colonyno});
    end
    
    
end
save(outfile, 'peaks');
end
