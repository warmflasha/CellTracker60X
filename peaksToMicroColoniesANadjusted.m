function [colonies,peaks]=peaksToMicroColoniesANadjusted(peaks)
% will have 7 columns to begin with, then the 8th column gets filled up
% here with dummy number(or can fill it in with the cytoplasmic area from IlastikPlusWatershed)
% and the 9th column gets filled in here with the colony size that each cell belongs to
% this way will run the tracker on these peaks and the data from colony
% grouping will be carried over to the cell structure ( the 5th column of cells.data)


global userParam;
userParam.colonygrouping = 120;

if ~exist('mm','var')
    mm=1;
end
for k=1:length(peaks)% fill in column 8 of peaks with dummy number
if ~isempty(peaks{k})
peaks{k}(:,8)=-1;
end
end
for ii=1:length(peaks)
    if ~isempty(peaks{ii})
                for k=1:size(peaks{ii},1)
            
               totcells = size(peaks{ii},1);
               %get number of columns from first non-empty image
                q=1; ncol=0;
                while ncol==0
                    ncol=size(peaks{ii},2);
                    q=q+1;
                end
                
                alldat=zeros(totcells,ncol+1);
                pts = [peaks{ii}(:,1) peaks{ii}(:,2)];
                alldat = peaks{ii};
                % pts should be the cells in the current image
                
                allinds=NewColoniesAW(pts);
                alldat = [alldat, allinds];
                ngroups = max(allinds);
                
                %Make colony structure for the single cell algorythm
                
                for jj=1:ngroups;
                    cellstouse=allinds==jj;
                    colonies(jj)=colony(alldat(cellstouse,:));%[1024 1344]
                end
            
        end
       
      
       
        % the code that fills up the 9th column of peaks with the size of the colony
        % that the cell belongs to
        peaks{ii}(:,9)=ones(size(peaks{ii},1),1); % initially assume that all the colonies are  of
        %size one (then the code fills out the ones that are not)
        ncel = [colonies.ncells]';%
        Ncol=size(ncel,1);% total colonies identified
        allcells=size(peaks{ii},1) ; % total number of cells within the peaks
        if Ncol == allcells
            peaks{ii}(:,9) = ones(size(peaks{ii},1),1);
        end
        
        %'allcells' are grouped into 'Ncol' colonies
        
        if allcells > Ncol
            [r,~]=find(ncel >1);
            [r2,~]=find(ncel == 1);
            onecel=allcells - size(r2,1); % number of one-cell colonies
            for j = 1:size(r,1)
                sz(j) = ncel(r(j)); % this is the size of the colony  the cell that is numbered
                %'r' in the peaks{ii} belongs to, so this cell and the (sz-1) cells after it get the label ncel(r) since all of them belong to the
                %colony of same size 'sz'
                peaks{ii}(r(j):(r(j)+sz(j)-1),9) = sz(j);
            end
            
            
        end
    end

end
%TO DO: add generic saving structure
 %save('OutfileFrame1','peaks','NucMasks','CytoMasks','colonies');% update
end