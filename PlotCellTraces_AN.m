
function [alldata] = PlotCellTraces_AN(dir, col,N)



[nums, filescyto]=folderFilesFromKeyword(dir,'Outfile1');


for j=1:length(nums)

matfile = filescyto(j).name;

pp=load(matfile,'peaks','cells'); % AN
% check the number of trajectories withn the frame and add the loop over
% same trajectories

%  trN = max(peaks{ii}(:,8));
%for k=1:trN
    
for ii=1:length(peaks)
          if ~isempty(peaks{ii})
       for k=1:length(peaks{ii},1)
           if peaks{ii}(k,9) == N && peaks{ii}(k,8)== 1 % here the number of the trajectory has to be
                alldata(q:(q+nlines(ii)-1))=colnew{ii}.data(:,col);
    end
end
   
alllines = sum(nlines);
alldata=zeros(alllines,nc);

q=1;
for ii=1:length(colnew)
    if ~isempty(colnew{ii})&& colnew{ii}.ncells == N
         if length(col)==1
            alldata(q:(q+nlines(ii)-1))=colnew{ii}.data(:,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
        else
            alldata(q:(q+nlines(ii)-1))=mean(colnew{ii}.data(:,col(1)))./mean(colnew{ii}.data(:,col(2)));% make sure the col numbers are correct
        end
        q=q+nlines(ii);
    end
    
end
% so need to get all N-cell colonies from all frames, belonging to the same
% trajectory


end




end

