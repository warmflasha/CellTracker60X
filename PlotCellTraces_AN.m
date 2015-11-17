
function [alldata] = PlotCellTraces_AN(dir, col,N)

%plot cell traces for specific colonies using the output from the runTracker
%EDS
% in that output last column of peaks = has the trajectory number that the
% cell belongs to
% in the colonies cell array( saved in the otfile) the last column is the
% number of colony within that frma that the cells belong to

% 1. read in the files in the directory, find the ones that have the
% 'outfile' srting in them
% 2. setup a loop to precess each file (it's peaks colonies and cells
% 3. group the output based on cell belonging to the same trajectory and
% the a colony of size N ( add the number of cells within the colony to the
% cells structure.
% plot the output for different colony sizes

[nums, filescyto]=folderFilesFromKeyword(dir,'Outfile');


for j=1:length(nums)

matfile = filesnuc(j).name;

pp=load(matfile,'peaks','colonies','cells'); % AN
% check the number of trajectories withn the frame and add the loop over
% same trajectories
% 
for k=1:length(peaks)
    if ~isempty(peaks{k})
colnew{k}.data = [colonies{k}.data peaks{k}(:,end)]; % last column of peaks should have the trajectory number
colnew{k}.ncells = [colonies{k}.ncells];    
   
    end
end
% after this loop should have the colnew cell array with each cell having
% all necesary data in the last two columns(trajectory number and number of
% cells within ceolony)
% each matfile name is overwritten

nc = N;% size(colnew{1}.data(end),2)number of cells within the colony 
alldata = zeros(length(colnew),nc);
nlines=zeros(length(colnew),nc);
    
for ii=1:length(colnew)
    if ~isempty(colnew{ii})&& colnew{ii}.ncells == N
                nlines(ii)=size(colnew{ii},1);
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

% pp=load(matfile,'peaks','NucMasks','colonies'); % AN
% peaks=pp.peaks;
% colonies = pp.colonies;
% nc = size(colonies,2);% number of colonies found in the image NEED TO GET THE MAX NUMBER OF COLONIES IN THE IMAGE
% alldata = zeros(length(colonies),nc);
% nlines=zeros(length(colonies),nc);
%     
% for ii=1:length(colonies)
%     if ~isempty(colonies{ii})&& colonies{ii}.ncells == N
%                 nlines(ii)=size(colonies{ii},1);
%     end
% end
%    
% alllines = sum(nlines);
% alldata=zeros(alllines(1),nc);
% 
% q=1;
% for ii=1:length(colonies)
%     if ~isempty(colonies{ii})
%          if length(col)==1
%             alldata(q:(q+nlines(ii)-1))=colonies{ii}.data(:,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
%         else
%             alldata(q:(q+nlines(ii)-1))=mean(colonies{ii}.data(:,col(1)))./mean(colonies{ii}.data(:,col(2)));
%         end
%         q=q+nlines(ii);
%     end
%     
% end


end

