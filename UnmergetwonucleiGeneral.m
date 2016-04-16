function [MaskFin2,maskfin] = UnmergetwonucleiGeneral(mask3)% input the nuclear mask
% the mask to input should be already thresholded binary from probability mask of ilastik
% function checks if the nuclear mask has two cells merged
% if yes, then splits them and returns a new mask where all the two merged
% nuclei are separated
% checks every nucleis in the mask, not limited to one object;
% if there are no merged objects then returns the same mask as input
clear MaskFin2
clear maskfin
clear masktmp

%
% get the image with only the merged object
global userParam;
mask3new = bwareafilt(mask3,[300 20000]);  %  userParam.areanuclow
stats = bwconncomp(mask3new);

nn = (stats.NumObjects);

if nn == 0
    MaskFin2 = mask3;
    maskfin = mask3;
    return
end

masktmp = cell(1, nn);

for ii=1:nn
    masktmp{ii} = zeros(1024,1024);
     MaskFin2{ii} = masktmp{ii};
    masktmp{ii}(stats.PixelIdxList{ii}) = 1;
    
end
for ii=1:nn
    [a,b,rx,ry,ch_x,ch_y,extrafilt,data_ch,data_mc,data_c2] = testmergenuclei(masktmp{ii});
    % disp(a);disp(b);
    if (a<60 || b<60)  || (a == 0) || (b == 0)             % refine this condition (not to split the nuc)
        MaskFin2{ii} = masktmp{ii};
        disp('no split');
        continue
    
    else
        
        [v] = checkboundaries(extrafilt,data_ch,data_c2);
        if v == 1
            MaskFin2{ii} = masktmp{ii};
            continue
        end
        [toelim2] = maxalongX(extrafilt,data_ch,data_mc,data_c2);
        [toelim2_y] = maxalongY(extrafilt,data_ch,data_mc,data_c2);
        
        
        if isempty(toelim2_y)&&isempty(toelim2)
            
            MaskFin2{ii} = masktmp{ii} ;
            disp('no extra');
             continue%
        end
        if isempty(toelim2)
            toelimfin = toelim2_y;
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',1));
            MaskFin2{ii} = masktmp{ii}&~II ;
            disp('one empty');
            didsplit = bwconncomp(MaskFin2{ii});
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area < 2000 || stats(2).Area < 2000)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc');
            end
            %continue
        end
        
        if isempty(toelim2_y)
            toelimfin = toelim2;
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',1));
            MaskFin2{ii} = masktmp{ii}&~II ; %
            disp('other empty');
            didsplit = bwconncomp(MaskFin2{ii});
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area < 2000 || stats(2).Area < 2000)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc');
            end
            %continue
        end
        
        if size(toelim2_y,1)>=size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',2)); % 'disk',4
            MaskFin2{ii} = masktmp{ii}&~II ;
            disp('split fine');
            didsplit = bwconncomp(MaskFin2{ii});
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area < 2000 || stats(2).Area < 2000)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc');
            end
            
        end
        if size(toelim2_y,1)>size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',2)); % 'disk',4
            MaskFin2{ii} = masktmp{ii}&~II ;
            didsplit = bwconncomp(MaskFin2{ii});
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==1     %
                
                toelimfin = toelim2_y;
                I = zeros(1024,1024);                                    % create an image with only that element
                linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
                I(linearInd)=1;
                II = imdilate(I,strel('disk',1)); % 'disk',4             % dilate a little in order to create a merged line
                
                %remove those pixels from the merged object
                MaskFin2{ii} = masktmp{ii}&~II; %
                disp('use toelim2_y');
            end
            %here check the area of the objects(in case split one nuc) 
            if didsplit.NumObjects ==2 && (stats(1).Area < 2000 || stats(2).Area < 2000)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc');
            end
            
        end
        if size(toelim2_y,1)<size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2_y;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',5)); % 'disk',4
            MaskFin2{ii} = masktmp{ii}&~II ;
            disp('split fine 2');
            didsplit = bwconncomp(MaskFin2{ii});
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area < 2000 || stats(2).Area < 2000)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc');
            end
            % continue
        end
        
    end
end

maskfin = zeros(1024,1024);
for ii=1:size(MaskFin2,2)
    I = MaskFin2{ii};
    stats = regionprops(I,'PixelIdxList');
    t = cat(1,stats.PixelIdxList);
    maskfin(t) = 1;
end
maskfin = im2bw(maskfin);
maskfin = bwareafilt(maskfin,[900 20000]);

end