function [colonies, peaks]=peaksToMicroColoniesAN(matfile,img)

%pp=load(matfile,'peaks','acoords','imgfiles','dims','userParam'); % AN
pp=load(matfile,'peaks','NucMasks','dims'); % AN

peaks=pp.peaks;
dims = size(peaks);
%ac=pp.acoords;
%dims=pp.dims;
%param = pp.userParam;
 
% if any(dims > 1)
% peaks=removeDuplicateCells(peaks,ac);
% end
if ~exist('mm','var')
    mm=1;
end

% k1=num2cell(ones(1,length(peaks)));
% lens=cellfun(@size,peaks,k1);
% totcells=sum(lens);
totcells = size(peaks{img},1);

%get number of columns from first non-empty image
q=1; ncol=0;
while ncol==0
    ncol=size(peaks{img},2);
    q=q+1;
end

alldat=zeros(totcells,ncol+1);

% q=1;
% for ii=1:length(peaks)
%    % if ~isempty(peaks{ii})
%         currdat=peaks{ii};
%         toadd=[ac(ii).absinds(2) ac(ii).absinds(1)];
%         currdat(:,1:2)=bsxfun(@plus,currdat(:,1:2),toadd);
%         alldat(q:(q+lens(ii)-1),:)=[currdat ii*ones(lens(ii),1)];
%         q=q+lens(ii);
%     end
% end
% pts=alldat(:,1:2);
    pts = [peaks{img}(:,1) peaks{img}(:,2)];
    alldat = peaks{img};
% pts should be the cells in the current image
   
    allinds=NewColoniesAW(pts);
    alldat = [alldat, allinds];
    ngroups = max(allinds);
    
    %Make colony structure for the single cell algorythm
    for ii=1:ngroups;
        cellstouse=allinds==ii;
        colonies(ii)=colony(alldat(cellstouse,:));%[1024 1344]
    end
    
    %put data back into peaks
    for ii=1:length(peaks)
        cellstouse=alldat(:,end-1)==ii;
        peaks{ii}=[peaks{ii} alldat(cellstouse,end-1:end)];
    end
    
    %peaks{img}=[peaks{img} alldat(cellstouse,end-1:end)];
    
end

