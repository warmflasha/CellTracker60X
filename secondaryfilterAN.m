function [zrange,smasks] = secondaryfilterAN(pmasks, minobjzstart, minsolid)

% zstart:first z slice with at least 4 objects; objects tracked from this slice.

for z = 1:size(pmasks,3)
    objects = regionprops(pmasks(:,:,z),'Centroid');
    if(size(objects,1)>=minobjzstart)
        zstart = z;
        break;
    end
end

zend = size(pmasks,3);

if (~exist('zstart'))
    smasks = [];
    zrange =0;
    
else
    cellsp = 1;
    
    zrange = [zstart:zend];
    
    smasks = false(size(pmasks));
    
%>>>>>>> upstream/master
    for z= zstart:zend
        
        tmp = pmasks(:,:,z);
        CC = bwconncomp(tmp);
        stats = regionprops(CC, 'solidity');
        
        m = 1;
        for solid = minsolid
            
            bad1{m} = find([stats.Solidity] < solid);
            tmp1{m}{z} = tmp;
            
            for i = 1:numel(bad1{m})
                tmp1{m}{z}(CC.PixelIdxList{bad1{m}(i)}) = 0;
            end
        tmp1{m}{z} = tmp1{m}{z};       
            
             %tmp1{m}{z} = bwareafilt(tmp1{m}{z}, [sizefilter sizefilter*50]);
            
            m = m+1;
        end
        
        
    end
% <<<<<<< HEAD
%   
%     
%     %%
%     % tmpn: modified clean masks 
% =======
%     
%     
%     %%
%     % tmpn: modified clean masks
% >>>>>>> upstream/master
    % adding unique objects from low sol to high sol
    % consz = no. of consecutive z slices in which unique objects are checked
    % overlap: only if the new unique object does not overlap with any previous object, it is added to the corresponding z slice.
    % zrange : z slices selected for further analyses
    
    consz = 2;
    overlap = 0;
    zrange = zstart:zend;

%     
%>>>>>>> upstream/master
    tmpn = soliditycombine(tmp1, consz, overlap, zrange);
    for z = zrange
        smasks(:,:,z) = tmpn{1}{z};
        
    end
    
end
% <<<<<<< HEAD
%     
% =======
end
%>>>>>>> upstream/master
