% plot distributions for live cell data
function [peaksnew1,peaksnew2] = histSignaling(fr_stim, resptime)


ff = dir('*_test*.mat');
q = 1;
peaksnew1 = [];
peaksnew2 = [];

for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
    
    % outfile = ('12_3D_20hr_test_xyz.mat');
    load(outfile,'colonies','peaks');
    tps = length(peaks);
  
    for j=1:fr_stim
        sizenew  = size(peaks{j},1);
        if ~isempty(peaks{j})
        peaksnew1(q:(q+sizenew-1),1) = peaks{j}(:,6)./peaks{j}(:,7);
        q = q+sizenew;
        end
    end
    for j=(fr_stim+resptime):tps
        sizenew  = size(peaks{j},1);
        if ~ isempty(peaks{j})
        peaksnew2(q:(q+sizenew-1),1) = peaks{j}(:,6)./peaks{j}(:,7);
        q = q+sizenew;
        end
    end
     
    
end
 peaksnew1(isinf(peaksnew1) == 1) = [];
 peaksnew1((peaksnew1) == 0) = [];
 peaksnew2(isinf(peaksnew2) == 1) = [];
 peaksnew2((peaksnew2) == 0) = [];
 histogram(peaksnew1,'Normalization','pdf');hold on
 histogram(peaksnew2,'Normalization','pdf'); 
 ylim([0 3]);
 xlim([0 3]);
 
end