
function [tr_new1,tr_new2] = histSignalingCol(timecolSZ, resptime,N,index,flag)
%input

%output

% get distribution plots for each colony size


fr_stim = 22;          %22 %38 %16
delta_t = 12; 
p = fr_stim*delta_t/60;
p2 = (timecolSZ)*delta_t/60;

%C = {'g','r','b','m'};
ff = dir('*_3Dsegm_febdata*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
clear traces 
clear traces_one
clear traces_two
clear traces_three

q = 1;
r = 1;
p = 1;

for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
  % outfile = ('12_3D_20hr_test_xyz.mat');
    load(outfile,'colonies','peaks');
    tps = length(peaks);
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    
    for j = 1:numcol
        colSZ = colonies(j).numOfCells(timecolSZ); % colony size determined at the time of stimulation
        if size(index,2) > 1   % look at the nuc to smad ratio
        traces{j} = colonies(j).NucSmadRatio(:); 
        end
        if size(index,2) == 1   % look at the nuc mean intensity in the rfp channel
         traces{j} = colonies(j).NucOnlyData(:);    
        end
        traces{j}(traces{j}==0) = NaN;
        traces(cellfun(@isempty,traces)==1)=[];
        if colSZ>0 && colSZ == N
            traces{j}(isnan(traces{j})==1) = 0;
            d =  size(traces{j},2);
%             if colSZ == 1
            traces_one{q}= traces{j};
%             end
%             if colSZ == 2
%             traces_two{r}= traces{j};
%             end
%             if colSZ == 3
%             traces_three{p}= traces{j};
%             end
        end
         q = q+1; 
         p = p+1;
         r = r+1;
         
    end
end

traces_one(cellfun(@isempty,traces_one)==1)=[];

d = size(traces_one,2);
clear replace
clear sm

for k=1:d
a =size(traces_one{k},1);
if a < 99;
sm(k) = a;
end
end
a = find(sm>0);
sm1 = nonzeros(sm);

for jj=1:size(nonzeros(sm),1)
    
    replace{jj} = zeros(99,1);
    replace{jj}(1:sm1(jj),1) = traces_one{a(jj)}(:,1);
    traces_one{a(jj)} = replace{jj};
end
clear traces_one_new
q = 1;
for k=1:d
s = size(traces_one{k},2);
traces_one_new(:,q:q+s-1) = traces_one{k}(:,:);
q = q+1;
end
q = 1;
h = 1;
K = size(traces_one_new,1)-fr_stim-resptime;
 for xx=1:size(traces_one_new,2)
     tr_new1(q:q+fr_stim-1,1) = traces_one_new(1:fr_stim,xx);
     tr_new2(h:h+K-1,1) = traces_one_new((fr_stim+resptime+1):end,xx);
     q = q+fr_stim;
     h = h+K;
 end
tr_new1(tr_new1==0)=[];
tr_new1(isinf(tr_new1)==1) = [];

tr_new2(tr_new2==0)=[];
tr_new2(isinf(tr_new2)==1) = [];
if flag == 1 && size(index,2)>1
histogram(tr_new1,'Normalization','pdf');hold on
histogram(tr_new2,'Normalization','pdf');
ylabel('Frequency');
xlabel('mean Nuc/Cyto smad4');
title(['All microColonies of size ' num2str(N) ]);
end
clear tr_new_all
if flag == 1 && size(index,2)==1
    tr_new_all = cat(1,tr_new1,tr_new2);
    histogram(tr_new_all,'Normalization','pdf');hold on
    ylabel('Frequency');
    xlabel('mean intensity in RFP channel');
    title(['All microColonies of size ' num2str(N) ]);
end

end