 % plot cell tracesfor the specific outfile
 
 function plotcelltraces(outfile,trajmin)
 %colormap = prism;
 C = {'b','r','g','m','b','c','g','r','b','g','m'};
 load(outfile,'colonies');
 % tps = length(peaks);
 numcol = size(colonies,2); % how many colonies were grouped within the frame
 traces = cell(1,numcol);
 for j = 1:numcol
     
     traces{j} = colonies(j).NucSmadRatio;%colonies(j).NucSmadRatio(:)
%      colSZ1 = colonies.ncells_predicted;
%      colSZ2 = colonies.ncells_actual;
     traces{j}((traces{j} == 0)) = nan;
     for h = 1:size(traces{j},2)
         a = isfinite(traces{j}(:,h));
         dat = zeros(size(traces{1},1),1);
         dat(a == 1,1) = traces{j}(a==1,h);
         if length(nonzeros(dat))>trajmin
             disp(['filter trajectories below' num2str(trajmin)]);
             disp(['use' num2str(length(nonzeros(dat)))]);
             hold on;figure(j+1),plot(dat,'-o','color',C{j});hold on% 
             %text(traces{j}(:,h)+0.1,num2str(colSZ1(:)));
             ylim([0 2.5]);
             ylabel('mean Nuc/Cyto smad4  ');
             xlabel('frames');
         end
     end
 end
 end
                 