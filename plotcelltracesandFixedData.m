 % plot cell tracesfor the specific outfile
 
 function plotcelltracesandFixedData(outfile,trajmin)
 % get time dependent data for a given .mat file ( plot each colony in
 % different figure)
 %colSZ = 3;
 %p2 = (resptime+jumptime)*delta_t/60;
 %coloniestoanalyze = 3;
 %cmap = summer;
 C = {'r','b','g','m','b','c','g','r','b','g','m'};
 %outfile = ('16_40X_tbsht1Z.mat');%3Dsegm_febdata    3D_20hr_test_xyz   % 7 ,16, 15, 18 29 26 31  reanalyze
 load(outfile,'colonies2','peaks');
 % tps = length(peaks);
 numcol = size(colonies2,2); % how many colonies were grouped within the frame
 traces = cell(1,numcol);
 l = 100;
 N = 0;
 for j = 1:numcol
     
     traces{j} = colonies2(j).NucSmadRatio;%colonies(j).NucSmadRatio(:)
     %colSZ = colonies(j).numOfCells(timecolSZ);
     traces{j}((traces{j} == 0)) = nan;
     for h = 1:size(traces{j},2)
         
         if length(traces{j}(isnan(traces{j}(:,h))==0))>trajmin
             figure(j), plot(traces{j}(:,h),'-*','color',C{j});hold on% cmap(j,:) 'r' traces
             
             ylim([0 2.5]);
             ylabel('mean Nuc/Cyto smad4  ');
             xlabel('frames');
             if colonies2(j).cells(h).onframes(end) == (l-N)
                 figure(j), hold on
                 text(colonies2(j).cells(h).onframes(end),traces{j}(end,h),num2str(colonies2(j).cells(h).fluorData(1,end)),'color','r','fontsize',15);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
             end
         end
     end
 end