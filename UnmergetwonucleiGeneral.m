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
clear trans
%
% get the image with only the merged object
global userParam;
 
mask3new = bwareafilt(mask3,[userParam.areanuclow userParam.areanuclow*10]);  %  userParam.areanuclow
stats = bwconncomp(mask3new);
stats2 = regionprops(mask3new,'Area');
mm = [stats2.Area];
maxel = mm(mm==max(mm));
%userParam.minnucfragment =2000;%  2500set this value to useless , so that the cutting is not done on whole  cells( move this choice to before cutting)
% decide whether to cut based on the values of a,b, below
userParam.areanuclow_unmerge = mean(mm) ;

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
     disp(a);disp(b);
    if (a < userParam.tocut || b < userParam.tocut)  || (a == 0) || (b == 0)             % 245 refine this condition (not to split the nuc)
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
        
        
        if isempty(toelim2_y) == 1 && isempty(toelim2)==1
            
            MaskFin2{ii} = masktmp{ii} ;
            disp('no extra');
             continue%
        end
        if isempty(toelim2) == 1
            toelimfin = toelim2_y;
            I = zeros(1024,1024);                                    % create an image with only that element
            
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',userParam.linedil));
            trans = masktmp{ii}&~II;
            didsplit = bwconncomp(trans);
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area <= userParam.minnucfragment || stats(2).Area <= userParam.minnucfragment)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc1');
            else
                MaskFin2{ii} = masktmp{ii}&~II ;
                disp('one empty');
            end
            %continue
        end
        
        if isempty(toelim2_y) == 1
            toelimfin = toelim2;
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',userParam.linedil));
            trans = masktmp{ii}&~II ; 
            didsplit = bwconncomp(trans);
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area <= userParam.minnucfragment || stats(2).Area <= userParam.minnucfragment)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc2');
            else
                MaskFin2{ii} = masktmp{ii}&~II ; %
            disp('other empty');
            end
            %continue
        end
        
        if size(toelim2_y,1)>=size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',userParam.linedil)); % 'disk',4
            trans= masktmp{ii}&~II;
            didsplit = bwconncomp(trans);
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area <= userParam.minnucfragment) || didsplit.NumObjects ==2 && (stats(2).Area <= userParam.minnucfragment)% 2300
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc3');
            else
                MaskFin2{ii} = masktmp{ii}&~II ;
                disp('split fine');
            end
            
        end
        if size(toelim2_y,1)>size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',userParam.linedil)); % 'disk',4
            trans = masktmp{ii}&~II ;
            didsplit = bwconncomp(trans);
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==1     %
                
                toelimfin = toelim2_y;
                I = zeros(1024,1024);                                    % create an image with only that element
                linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
                I(linearInd)=1;
                II = imdilate(I,strel('disk',userParam.linedil)); % 'disk',4             % dilate a little in order to create a merged line
                
                %remove those pixels from the merged object
                MaskFin2{ii} = masktmp{ii}&~II; %
                disp('use toelim2_y');
            end
            %here check the area of the objects(in case split one nuc) 
            if didsplit.NumObjects ==2 && (stats(1).Area <= userParam.minnucfragment || stats(2).Area <= userParam.minnucfragment)
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc4');
            else
                MaskFin2{ii} = masktmp{ii}&~II ;
                disp('use toelim2');
            end
            
        end
        if size(toelim2_y,1)<size(toelim2,1) && ~isempty(toelim2) && ~isempty(toelim2_y)
            toelimfin = toelim2_y;
            
            I = zeros(1024,1024);                                    % create an image with only that element
            linearInd = sub2ind(size(I), toelimfin(:,2), toelimfin(:,1));
            I(linearInd)=1;
            II = imdilate(I,strel('disk',userParam.linedil)); % 'disk',4
            trans = masktmp{ii}&~II ;
            
            didsplit = bwconncomp(trans);
            stats = regionprops(didsplit,'Area','Centroid');
            if didsplit.NumObjects ==2 && (stats(1).Area <= userParam.minnucfragment) || didsplit.NumObjects ==2 && (stats(2).Area <= userParam.minnucfragment)% 2300
                MaskFin2{ii} = masktmp{ii} ;
                disp('no split, single nuc5');
            else
                MaskFin2{ii} = masktmp{ii}&~II ;
            disp('split fine 2');
            end
            % continue
        end
        
    end
end

maskfin = zeros(1024,1024);

for ii=1:size(MaskFin2,2)
    J = MaskFin2{ii};
    stats = regionprops(J,'PixelIdxList');
    t = cat(1,stats.PixelIdxList);
    maskfin(t) = 1;
end
maskfin = im2bw(maskfin);
maskfin = bwareafilt(maskfin,[userParam.areanuclow2 userParam.areanuclow2*10]);

end