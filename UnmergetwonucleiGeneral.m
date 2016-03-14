function [MaskFin2,maskfin] = UnmergetwonucleiGeneral(mask3)% input the nuclear mask
% the mask to input should be already thresholded binary from probability mask of ilastik
% function checks if the nuclear mask has two cells merged
% if yes, then splits them and returns a new mask where all the two merged
% nuclei are separated
% checks every nucleis in the mask, not limited to one object;
% if there are no merged objects then returns the same mask as input
clear MaskFin2
clear maskfin
%
% get the image with only the merged object
global userParam;
mask3new = bwareafilt(mask3,[600 20000]);  %  userParam.areanuclow
stats = bwconncomp(mask3new);

nn = (stats.NumObjects);
masktmp = cell(1, nn);

for ii=1:nn
    masktmp{ii} = zeros(1024,1024);
    MaskFin2{ii} = masktmp{ii};
    masktmp{ii}(stats.PixelIdxList{ii}) = 1;
    
end
for ii=1:nn
    [a,b,extrafilt,data_ch,data_mc,data_c2] = testmergenuclei(masktmp{ii});
    if (a<30 || b<30) 
        MaskFin2{ii} = masktmp{ii};
    else

[v] = checkboundaries(extrafilt,data_ch,data_c2);
if v == 1
  MaskFin2{ii} = masktmp{ii};
  continue
end
        [toelim2] = maxalongX(extrafilt,data_ch,data_mc,data_c2);
        [toelim2_y] = maxalongY(extrafilt,data_ch,data_mc,data_c2);
        
        
        % if isempty(hull_bnd_only1)||isempty(hull_bnd_only2)||isempty(obj_bnd_only1)||isempty(obj_bnd_only2)
        %   MaskFin2{ii} = masktmp{ii};
        %   continue
        % end
        
        % Choose which line has the least intersecting points with the merged
        % object = this is the line where the cut wll be made
        % if ((objrow2(rX2(1)))~= 0 && objrow(rX(1))~=0)==0
        %     toelimfin = toelim2_y;
        %     I = zeros(1024,1024);
        %     linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
        %     I(linearInd)=1;
        %     II = imdilate(I,strel('disk',4)); % 'disk',4
        %
        %     MaskFin2{ii} = masktmp{ii}&~II ;  %
        %     %MaskFin2 = MaskFin + mask3old;
        % else
        if isempty(toelim2_y)&&isempty(toelim2)
            
             MaskFin2{ii} = masktmp{ii} ; %
         end 
        if isempty(toelim2)
            toelimfin = toelim2_y;
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',4));
             MaskFin2{ii} = masktmp{ii}&~II ; %
        end

         if isempty(toelim2_y)
            toelimfin = toelim2;
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',4));
            MaskFin2{ii} = masktmp{ii}&~II ; %
         end
        
         if size(toelim2_y,1)>=size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',4)); % 'disk',4
            
            didsplit = bwconncomp(II);
            if didsplit.NumObjects ==1 
                toelimfin = toelim2_y;
            end
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',4)); % 'disk',4             % dilate a little in order to create a merged line
                  
            
                                   % remove those pixels from the merged object
            MaskFin2{ii} = masktmp{ii}&~II ; %  
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
end