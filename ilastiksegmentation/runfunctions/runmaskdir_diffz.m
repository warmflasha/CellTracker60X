function runmaskdir_diffz(segfiledir, rawfiledir, mrnafilepath, paramfile, outfiledir, fluorpdir)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% mrnafilepath : directory path of spatzcells output.
% nzslices: no. of z slices
% outfile: final output file path

global userparam
eval(paramfile);
mkdir (outfiledir);


ff = readFISHdir(rawfiledir, userparam.nsamples);
nuc_ch = userparam.nucchannel;

for samplenum = 1:userparam.nsamples
    
    npositions = ff.positions(samplenum);
    
    clear peaks1 peaks;
    %for pos = 12;
    for pos = 1:npositions
        nzslices = ff.zslices{samplenum}(pos);
        
        
        [pnuc, inuc] = readmaskfilesnew(segfiledir, rawfiledir, samplenum, pos-1,  nzslices, nuc_ch);
        %%
        
        pmasks = primaryfilter(pnuc,userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
        %%
        
        [zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
        
        if (zrange)
            %%
            zmatch = nzslices-1;
            [PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userparam.matchdistance, zrange, zmatch);
            %%
            
            
            [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);
            
            %%
            
            framenum = sum(ff.positions(1:samplenum-1)) + pos;
            peaks1{pos} = mrnapercells(nucleilist, stats, mrnafilepath, framenum, zrange, userparam.channels, userparam.cmcenter);
            %%
            peaks{pos} = peakscelltrackerformat(peaks1{pos});
            
        else
            peaks{pos} = [];
        end
        
        %%
        if(exist(fluorpdir) && samplenum ~= userparam.negativecontrol)
            [dirinfo2, start2] = readdirectory(fluorpdir);
            
            if samplenum < userparam.negativecontrol
                fimageno = start2+(samplenum-1)*nzslices;
            else
                fimageno = start2+(samplenum-2)*nzslices;
            end
            peaks{samplenum} = masksin3d(CC, nucleilist, nzslices, fimageno, fluorpdir, dirinfo2, inuc, peaks{samplenum});
        end
        
        
    end
    
    outfile = strcat(outfiledir, '/', sprintf('sample%02dout.mat', samplenum));
    save(outfile, 'peaks');
end
