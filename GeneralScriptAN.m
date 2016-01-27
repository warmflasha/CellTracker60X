%%
%ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
ilastikfile = ('/Users/warmflashlab/Desktop/Dec31setIlastikMasks_headless_DiffW0/NucMaskDiffDec31set_tg3.h5');
ilastikfile2 = ('/Users/warmflashlab/Desktop/Dec31setIlastikMasks_headless_DiffW1/CytoMaskDiffDec31set_tg3.h5');
nuc = h5read(ilastikfile,'/exported_data');
cyto = h5read(ilastikfile2,'/exported_data');
  k =1;
    nuc = nuc(2,:,:,k);% for probabilities exported
    nuc = squeeze(nuc);
    mask1 = nuc;
    
    cyto = cyto(2,:,:,k);% for probabilities exported
    cyto = squeeze(cyto);
    mask2 = cyto;
    
    figure(2), subplot(1,2,1),imshow(mask1);
   figure(2), subplot(1,2,2),imshow(mask2);
    Lnuc = im2bw(mask1,0.7);
   Lcyto = im2bw(mask2,0.7);
   figure, imshow(Lcyto'&~Lnuc');
%

%%
ff = dir('*jan8set*.mat');
k = 1;
outfile = ff(k).name;
load(outfile);

%%
N = 34;
n = uncompressBinaryImg(imgfiles(N).compressNucMask);
nc = uncompressBinaryImg(imgfilescyto(N).compressNucMask);
%close all
figure,subplot(1,2,1),imshow(n);

hold on
subplot(1,2,2),imshow(nc);
%%
fr_stim = 16;
ff = dir('*jan8set*.mat');
for k=1:size(ff,1)
    outfile = ff(k).name;
    
runTracker(outfile,'newTrackParamAN');
global userParam;
userParam.colonygrouping = 120;
% look at colonies around the stimulation frame (window of couple hours)?
cellsToDynColonies(outfile);
end

%%
% colonies(2).numOfCells(1)
ff = dir('*Outfile*.mat');%jan8set % Pluri_42hrs % 22hr set: k = 3,4,6,16,22 traces k = 10 bad  %Outfile
k =6;
outfile = ff(k).name;
load(outfile);

fr_stim = [];%22
fldat = [2 3];
delta_t = 5; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime = 36;% 36in frames ( converted to hours later)
GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);

%%
for j=1:size(ncells,1)
plot(ncells{j},'-*');hold on
ncells{j}(10)
end

%%
%OLD:  plot the stats (before, after, amplitude...)
fr_stim = 22;%22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime =50;% 15 50 36 in frames ( converted to hours later)
coloniestoanalyze = 3; % max size of the colonies that want to look at
for jj = 1:coloniestoanalyze
    
    colSZ = jj;
    
    %GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);
    %ff = dir('*10ngmlDifferentiated*.mat');
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs
    for k=1:size(ff,1)
        outfile = ff(k).name; %nms{k};
        load(ff(k).name);
        datafin = GetDynamicColonyStats(outfile,fr_stim,delta_t,flag,colSZ,resptime,coloniestoanalyze);
        
    end
end
if ~isempty(fr_stim)
figure(11)
xx = 0:0.1:2.5;
plot(xx,xx,'-k','linewidth',2);

end
if isempty(fr_stim) 
 figure(11) 
 xx = 0:1:coloniestoanalyze;
 yy = ones(1,coloniestoanalyze+1);
 plot(xx,yy,'--');
 figure(12)
 plot(xx,yy,'--');
end
%%
% NEW: get new traces
fr_stim = 1 ;%22 %38 %16
fldat = [2 3];
delta_t = 12; % 12 %15  in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;
colSZ = 2;
resptime =50;% 15 50 36 in frames ( converted to hours later)
jumptime = 4;% in frames
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = summer;
C = {'b','g','r','m'};
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
    %for %k=1:size(ff,1);
    k=29;

        outfile = ff(k).name; %nms{k};
        %cellsToDynColonies(outfile);
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        traces = cell(1,numcol);
        
        for j = 1:numcol
        
          
            %tpt =  (colonies(j).cells(k).onframes');
            traces{j} = colonies(j).NucSmadRatio(:);
                       
            figure(j), plot(traces{j},'-*','color',cmap(j,:));% cmap(j,:) 'r' traces
            ylim([0 2]);
            ylabel('mean Nuc/Cyto smad4  ');
            xlabel('frames');
         
         end
   % end 
%%
% NEW: get the before and after  plots

fr_stim = 22;% 16 22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;
bf_fin = [];
aft_fin = [];
window_fin = [];
resptime =170;% 15 50 36 in frames ( converted to hours later)
range = [26 5];
jumptime = 4;% in frames
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = parula;
C = {'b','g','r','m'};
q = 1;
w = 1;
s = 1;

for jj = 1:coloniestoanalyze
    
    colSZ = jj;
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
    for k=1:size(ff,1)
        
        outfile = ff(k).name ;
        %cellsToDynColonies(outfile);
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        for j = 1:numcol
            if ~isempty(fr_stim)
                
                if colSZ == colonies(j).numOfCells(fr_stim); % how many cells within colony at time fr_stim
                    
                    stats =  colonies(j).DynNucSmadRatio(tps,fr_stim,resptime,range,jumptime);%,resptime
                    bf = (stats(:,1));
                    aft =(stats(:,2));
                    window = (stats(:,3));
                    jump = (stats(:,4));
                    
                    b = find(isnan(bf));
                    bf(b)=0;
                    a = find(isnan(aft));
                    aft(a)=0;

                    w = find(isnan(window));
                    window(w)=0;
                   
                   currlengthbf = size(nonzeros(bf),1);
                   currlengthaft = size(nonzeros(aft),1);
                   currlengthgap = size(nonzeros(window),1);
                   
                    bf_fin{jj}((q:(q+currlengthbf)-1),1) = nonzeros(bf);%sum(bf);
                    aft_fin{jj}((w:(w+currlengthaft)-1),1) = nonzeros(aft);%sum(aft);
                    window_fin{jj}((s:(s+currlengthgap)-1),1) = nonzeros(window);%sum(window);
%                     bf_fin{jj}(q,2) = size(nonzeros(bf),1);
%                     aft_fin{jj}(q,2) = size(nonzeros(aft),1);
%                     window_fin{jj}(q,2) = size(nonzeros(window),1);
                    
                    q = q+currlengthbf;
                    w = w+currlengthaft;
                    s = s+currlengthgap;
                    
                    hold on,figure(4), plot(bf,aft,'*','color',C{colSZ},'markersize',15);
                    title(['mean Nuc/Cyto smad4 ' num2str(p2) 'hrs after bmp4']);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,1), plot(colSZ,bf,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,2), plot(colSZ,aft,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    hold on,figure(10), plot(colSZ,window,'*','color',C{colSZ},'markersize',15);%bf,window
                    title(['mean Nuc/Cyto smad4 between ' num2str((range(1)*delta_t)/60) 'and' num2str((range(2)*delta_t)/60) 'hours']);
                ylim([0 2]);
                xlim([0 2]);
                  hold on,figure(11), plot(colSZ,jump,'*','color',C{colSZ},'markersize',15);%amplitude of the actual jump during the jumptime
                    title(['Amplitude of the jump, bf to ' num2str(((fr_stim+jumptime+range(2))*delta_t)/60) 'hours']);
                ylim([0 2]);
                xlim([0  (coloniestoanalyze+1)]);

                
                end
                save('meansdiff.mat','bf_fin','aft_fin','window_fin');
            end
            
            if isempty(fr_stim)
                if colSZ == colonies(j).numOfCells(10);%
                    stats =  colonies(j).DynNucSmadRatio(tps,fr_stim,resptime,range,jumptime);%,resptime
                    bf = (stats(:,1));
                    window = (stats(:,3));
                    b = find(isnan(bf));
                    bf(b)=0;
                    w = find(isnan(window));
                    window(w)=0;
                    
                    currlengthbf = size(nonzeros(bf),1);
                    currlengthgap = size(nonzeros(window),1);
                     
                    bf_fin{jj}((q:(q+currlengthbf)-1),1) = nonzeros(bf);%sum(bf);
                    window_fin{jj}((s:(s+currlengthgap)-1),1) = nonzeros(window);%sum(window);
                   
                    q = q+currlengthbf;
                    w = w+currlengthgap;
                   
                    hold on,figure(4), plot(colSZ,bf,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    xx = 0:1:4;
                    yy = ones(1,5);
                    figure(4),plot(xx,yy,'-k','linewidth',2);
%                     hold on,figure(10), plot(colSZ,window,'*','color',C{colSZ},'markersize',15);
%                     title(['mean Nuc/Cyto smad4 between ' num2str((range(1)*delta_t)/60) 'and' num2str((range(2)*delta_t)/60) 'hours']);
%                     ylim([0 2]);
%                     xlim([0 2]);
                end
                save('meanspluri.mat','bf_fin');
            end
            
        end
        
    end
    
    
end
if ~isempty(fr_stim)
                xx = 0:1:4;
                yy = ones(1,5);
                hold on,figure(4),plot(xx,xx,'-k','linewidth',2);
                ylim([0 2]);
                xlim([0 2]);
                hold on,figure(10),plot(xx,xx,'-k','linewidth',2);
                ylim([0 2]);
                xlim([0 2]);
                ylabel(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
                xlabel('mean Nuc/Cyto smad4 Before stimulation');
                hold on,figure(5),subplot(1,2,2),plot(xx,yy,'-k','linewidth',2);
                ylabel(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
                ylim([0 2]);
end
if isempty(fr_stim)
                hold on,figure(5),subplot(1,2,1),plot(xx,yy,'-k','linewidth',2);
                 ylim([0 2]);
                xlim([0 2]);
end
%%
% plot the mean values
jumptime = 4;% in frames
fr_stim = 22;%22 %38 %16
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;
test = cell(1,20);
%resptime =50;% 15 50 36 in frames ( converted to hours later)
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = parula;
C = {'b','g','r','m'};
%test = cell(1,coloniestoanalyze);
for jj = 1:coloniestoanalyze
    colSZ = jj;
    if isempty(fr_stim)
        load('meanspluri.mat');
        before(jj) = mean(bf_fin{jj});
               
        errbf(jj) = std(nonzeros(bf_fin{jj}));
        
        figure(6), errorbar(jj,before(jj),errbf(jj),'*','markersize',15,'color',C{jj});hold on
        ylabel('Mean Nuc/Cyto smad4');
        ylim([0.5 1.7]);
        xlim([0 (coloniestoanalyze+1)]);
    end
    if ~isempty(fr_stim)
        load('meansdiff.mat');
       
        before(jj) = mean(nonzeros(bf_fin{jj})); %sum(bf_fin{jj}(:,1))/sum(bf_fin{jj}(:,2));
        after(jj) = mean(nonzeros(aft_fin{jj}));%sum(aft_fin{jj}(:,1))/sum(aft_fin{jj}(:,2));
            
        errbf(jj) = std(nonzeros(bf_fin{jj}));
        erraft(jj) = std(nonzeros(aft_fin{jj}));
   
        figure(7), errorbar(jj,before(jj),errbf(jj),'*','markersize',15,'color',C{jj});hold on;
        figure(7), errorbar(jj,after(jj),erraft(jj),'.','markersize',20,'color',C{jj});hold on;
         legend('before','after');
         title(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
        ylabel('Mean Nuc/Cyto smad4');
        ylim([0.3 1.9]);
        xlim([0 (coloniestoanalyze+1)]);
        figure(8), errorbar(jj,after(jj),erraft(jj),'.','markersize',20,'color',C{jj});hold on;
         title(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
        ylabel('Mean Nuc/Cyto smad4');
        ylim([0.3 1.9]);
        xlim([0 (coloniestoanalyze+1)]);
    end
end

      
      
      
