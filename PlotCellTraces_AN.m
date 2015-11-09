
function [alldata] = PlotCellTraces_AN(matfile, col,N)
% plot live cell analysis results for a given colony size
% col = which column of peaks you want to analyze
% N = the size of the colonies you want to look at
% matfile = outall file, separate for each position(but all time points are
% included)


pp=load(matfile,'peaks','NucMasks','colonies'); % AN
peaks=pp.peaks;

colonies = pp.colonies;

nc = size(colonies,2);% number of colonies found in the image NEED TO GET THE MAX NUMBER OF COLONIES IN THE IMAGE
%sz = size(colonies{1},1);% size of the colonies found in the image
%colors = colorcube(ncell);
alldata = zeros(length(colonies),nc);
%alldata2 = zeros(length(colonies),nc);
%colors = {'r','b'};
nlines=zeros(length(colonies),nc);
    for k=1:nc


for ii=1:length(colonies)
    if ~isempty(colonies{ii})&& colonies{ii}(k).ncells == N
                nlines(ii,k)=size(colonies{ii}(k),1);
    end
end
    end

alllines = sum(nlines);
alldata=zeros(alllines(1),nc);

for k=1:nc


q=1;

for ii=1:length(colonies)
    if ~isempty(colonies{ii})
        
        if length(col)==1
            alldata(q:(q+nlines(ii)-1),k)=colonies{ii}(k).data(:,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
        else
            alldata(q:(q+nlines(ii)-1),k)=mean(colonies{ii}(k).data(:,col(1)))./mean(colonies{ii}(k).data(:,col(2)));
%             
        end
        q=q+nlines(ii);
    end
    
end

end
end

