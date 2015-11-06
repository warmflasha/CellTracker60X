function [colonies]=peaksToMicroColoniesAN(peaks)

 %pp=load(matfile,'peaks','NucMasks'); % AN

%peaks=pp.peaks;
% colonies = cell(length(peaks));
if ~exist('mm','var')
    mm=1;
end

for k=1:length(peaks)
    
if ~isempty(peaks{k}) 
 

totcells = size(peaks{k},1);

%get number of columns from first non-empty image
q=1; ncol=0;
while ncol==0
    ncol=size(peaks{k},2);
    q=q+1;
end

alldat=zeros(totcells,ncol+1);
pts = [peaks{k}(:,1) peaks{k}(:,2)];
alldat = peaks{k};
% pts should be the cells in the current image

allinds=NewColoniesAW(pts);
alldat = [alldat, allinds];
ngroups = max(allinds);

%Make colony structure for the single cell algorythm
for ii=1:ngroups;
    cellstouse=allinds==ii;
    colonies{k}(ii)=colony(alldat(cellstouse,:));%[1024 1344]
end
    
 %save('OutfileFrame1','peaks','NucMasks','CytoMasks','colonies');
 
% put data back into peaks NEED TO ADJUST for 60X microcolonies

%     for ii=1:size(peaks{img},1)
%         cellstouse=alldat(:,end-1)==ii;
%         peaks{img}=[peaks{img} alldat(cellstouse,end-1:end)];
%     end
end
end
end


