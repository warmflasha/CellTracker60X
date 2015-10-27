% combined analysis of 60X images so far
% frame 0: two one-cell colonies
% frame 19: two 0ne-cell colonies and one 2-cell colony
% frame 10: three one-cell colonies 
% frame 1: two 2-cell colonies 
% frame 13: one cell dividing into 2 (frames 12-16 show division)

dir = '.';
load('Frame00_Analysis.mat');
fr0col1 = onecellcol1;
fr0col2 = onecellcol2;
load('Frame10_Analysis.mat');
fr10col1 = onecellcol1;
fr10col2 = onecellcol2;
fr10col3 = onecellcol3;
load('Frame19_Analysis.mat');
fr19col1 = onecellcol1;
fr19col2 = onecellcol2;
% load('Frame13_Analysis.mat');
% for i=1:11
% fr13col1{i} = nucmeanInt{i}./cytomeanInt{i};
% end
% for i=17:35
% fr13col2{i} = nucmeanInt{i}./cytomeanInt{i};
% end
for k=3:34
    plot(k,fr0col1{k},'r--*'),hold on
    plot(k,fr0col2{k},'g--*'),hold on
    plot(k,fr10col1{k},'b--*'),hold on
    plot(k,fr10col2{k},'y--*'),hold on
    plot(k,fr10col3{k},'m--*'),hold on
    plot(k,fr19col1{k},'k--*'),hold on
    plot(k,fr19col2{k},'b--*');hold on
ylim([0.1 1.9]);
end
xlabel('time, frames');
 ylabel('nuc/cyto mean for cells in the frame');
 title('GFPsmad4RFPh2b cells, 10ng/ml bmp4 added after frame 16, ~ 9 hours imaging time');

% hold on
% plot(fr13col1,'*y');
% hold on
% plot(fr13col2,'*y');