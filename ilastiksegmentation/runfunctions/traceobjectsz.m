  %%
  function [PILsn, PILsSourcen, masterCCn, stats, nuclein1, zrange,CC] = traceobjectsz(smasks, matchdistance, zrange, zmatch)
    
    % tracking objects
    %very simple matching
    
    clear nuclein1 PILsn PILsSourceN overlayn
    n = zmatch;
    
    CC = {};
    stats = {};
    zstart = zrange(1);
    zend = zrange(end);
    
    
    
    % recalculating zstart
    
    for z = zrange(1):zrange(end)
        objects = bwconncomp(smasks(:,:,z));
        if(objects.NumObjects > 3);
            zstart = z+1;
            break;
        end
    end
    
    
    
    for z = zrange(1):zrange(end)               %%%%%%%%%
        
        CC{z} = bwconncomp(smasks(:,:,z));
        stats{z} = regionprops(CC{z}, 'Centroid', 'Area');
    end
    
    
    
    nSlices = zrange(end)-zrange(1)+1;
    nuclein1 = nan([CC{zrange(1)}.NumObjects nSlices]);
    dmax = matchdistance;
    
    % master PixelIdxList
    PILsn = CC{zrange(1)}.PixelIdxList;
    
    % keep track of which slice added an object to the PILs
    PILsSourcen = ones([numel(PILsn) 1])*zrange(1);
    
    for i = 1:CC{zrange(1)}.NumObjects
        nuclein1(i,1) = i;
    end
  %%  
    % Irrespective of the interval check size,  the logic for matching cells of two frames
    % will be the same.
    
    for z = zrange(1)+1
        clear CM1 CM2 D
        
        CM1 = cat(1,stats{z-1}.Centroid);
        CM2 = cat(1,stats{z}.Centroid);
        
        D = distmat(CM1, CM2);
        
        [row,col] = find(D < dmax);
        
        % which of the new slice are matched to the previous
        matched = sum(D < dmax, 1)';
        
        % if matched add to current track
        for j = 1:sum(matched)
            trackIdx = find(nuclein1(:,z-zrange(1)) == row(j));
            nuclein1(trackIdx,z-zrange(1)+1) = col(j);
            PILsn{trackIdx} = cat(1,PILsn{trackIdx}(:), CC{z}.PixelIdxList{col(j)}(:));
        end
        
        % if not matched make new track
        for j = 1:CC{z}.NumObjects
            if ~matched(j)
                newtrack = nan([1 zrange(end)-zrange(1)+1]);
                newtrack(z-zrange(1)+1) = j;
                nuclein1 = cat(1,nuclein1, newtrack);
                PILsn = [PILsn, CC{z}.PixelIdxList(j)];
                PILsSourcen = cat(1, PILsSourcen, z);
            end
        end
    end
   %%
    for z = zrange(1)+2:zrange(end)
        
        if(z<zrange(1)+n)
            frchksize = z-zrange(1);
            framestart = 1;
            
        else
            frchksize = n;
            framestart = z-n-zrange(1)+1;
        end
        
        
        clear CM dm row col mm c0 c00 c1;
        
        for i = 0:frchksize
            CM{i+1} = cat(1,stats{z-frchksize+i}.Centroid);
        end
        
        for i = 1:frchksize
            dm{i} = distmat(CM{i}, CM{frchksize+1});
        end
        
        for i = 1:frchksize
            [row{i}, col{i}] = find(dm{i} < dmax);
            mm{i} = sum(dm{i} < dmax, 1)';
        end
        
        frn = framestart;
        frnlim = z-zrange(1)+1; 
        mfrn = 1;
        

        %%
        while frn < frnlim 
            
            if(mfrn ==1)
                c1 = find(mm{mfrn}(:,1)==1); % matched with frame1
                c0 = find(mm{mfrn}(:,1)==0); % did not match with frame1
            end
            
            %%
            if(exist('c1', 'var')) 
                for i1 = c1'
                    clear olabel
                    olabel = find(col{mfrn}(:,1) == i1);
                    zsearch = frn;
                    clear trackIdx
                    trackIdx = find(nuclein1(:,zsearch) == row{mfrn}(olabel,1));
                    nuclein1(trackIdx, z-zrange(1)+1) = i1;
                    PILsn{trackIdx} = cat(1,PILsn{trackIdx}(:), CC{z}.PixelIdxList{i1}(:));
                end
            end
            clear c1;
            %%
            if(exist('c0', 'var') && (~isempty(c0)))
                if(frn < z-zstart)
                    clear c1 c00;
                    m1 = 1;
                    m2 = 1;
                    for i1 = c0'
                        if (mm{mfrn+1}(i1,1) == 1)
                            c1(m1,1) = i1;
                            m1 = m1+1;
                        else
                            c00(m2,1) = i1;
                            m2 = m2+1;
                        end
                    end
                    clear c0;
                    if(exist('c00', 'var'))
                        c0 = c00;
                    end
                end
            end
            %% assign a new track for elements that did not not match with any of the previous frames.
            if(frn == z-zstart)&& exist('c0', 'var')&& (~isempty(c0))
                for i1 = c0'
                    clear trackIdx;
                    trackIdx = size(nuclein1, 1)+1;
                    nuclein1(trackIdx, 1:nSlices) = NaN;
                    nuclein1(trackIdx, z-zrange(1)+1) = i1;
                    PILsn = [PILsn, CC{z}.PixelIdxList(i1)];
                    PILsSourcen = cat(1, PILsSourcen, z);
                end
                
            end
            
            frn = frn+1;
            mfrn = mfrn+1;
        end
    end
    
   %% 
    masterCCn = CC{zrange(1)};
    masterCCn.NumObjects = numel(PILsn);
    masterCCn.PixelIdxList = PILsn;
    zrange = [zrange(1):zrange(end)];