% code to transform the paeks nto the coordinates of the montage nd then
% find the colonies
%%
% load outfile
% load the last peaks and make an image with peaks(:,1:2) as the centroids of objects/ or even colonies; make a directory of such images named as Andor
% then make a montage out of those peaks and transform it as above)
 
% direc2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/lastTlivecellset/');
% [ac, fi]=alignManyPanelsAndorZstackMontage(direc2,[10 4],0,1);%
% ff = dir('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
% [ac2, fi2]=alignManyPanelsAndorZstackMontage('good_Backtocoord40xBF2_20160709_53148 PM',[7 11],0,1);
 
% from DATAStorage
 
% direc2 = ('/Volumes/data2/Anastasiia/totestClonyGrouping/lastTlivecellset');
% [ac, fi]=alignManyPanelsAndorZstackMontage(direc2,[10 4],0,1);%
% ff = dir('/Volumes/data2/Anastasiia/totestClonyGrouping/FinalOutfiles_completeTraces');
% [ac2, fi2]=alignManyPanelsAndorZstackMontage('good_Backtocoord40xBF2_20160709_53148 PM',[7 11],0,1);
 
load('/Volumes/data2/Anastasiia/totestClonyGrouping/align.mat');
 
%direc3 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/PeaksIntoMontage_lastT/');
%imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/nuc_raw');%rawimages_nuc%
t = imread('/Volumes/data2/Anastasiia/totestClonyGrouping/test4mid.tif');
%t = imread('test4mid.tif');
imshow(t(:,:,2),[]); % 2 - green chanel, 1 - red(4x10)
imcontrast
%%
pl = 2;
positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 0;
%ff2 = readAndorDirectory(imagedir1);
chan = [1];
chanal = 1;
timegroup = [];
l = 100;
%%
xyall = [];
for k=1:size(positions,2) % same as the loop over image numbers
    
    outfile = [ num2str(positions(k)) strdir ];
    
    load(outfile,'peaks');
    
    if ~isempty(peaks{l-N})       
        disp([int2str(k) ': ' int2str(size(peaks{l-N},1)) ' cells.']);
            
    xynow = bsxfun(@plus,peaks{l-N}(:,1:2),ac(k).absinds([2 1]));
    xyall = [xyall; xynow];
    end
end
%%
load('/Volumes/data2/tform.mat');
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [934-815, 1660-1330];%empirically determined from image
%toTranslate2 = 4*toTranslate-[377 424]; %correction empirically determined. [200 300]
toTranslate2 = [-140, 40];
%%
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);

xyall_rot = bsxfun(@minus,xyall, mean(xyall)); % mean subtract;
xyall_rot = xyall_rot*rotmat; %rotate about center
xyall_rot = bsxfun(@plus,xyall_rot,midpoint([2 1])); %add back

xyall_rot_trans = bsxfun(@plus,xyall_rot,toTranslate2);
%%
figure; imshow(t(:,:,1),[0 5e3]); hold on;
plot(xyall_rot(:,1),xyall_rot(:,2),'r.');
plot(xyall(:,1),xyall(:,2),'g.');
plot(xyall_rot_trans(:,1),xyall_rot_trans(:,2),'c.');

%%
% directory with segmented files
fixedsegm = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-09-20-nucmasks_backtocoord/Segmented_backtoCoordTile1';
strdirfix = '_testBacktocoord';
ff = dir(fixedsegm);
xyallfix = [];
positions2 = (0:76); % saved outfiles starting from 0


for k=1:size(positions2,2) % same as the loop over image numbers
    
    outfile = [ num2str(positions2(k)) strdirfix ];
    
    load(outfile,'peaks');
       
    if ~isempty(peaks{1})   % only one time point here    
        disp([int2str(k) ': ' int2str(size(peaks{1},1)) ' cells.']);
            
    xynowfix = bsxfun(@plus,peaks{1}(:,1:2),ac2(k).absinds([2 1]));
    xyallfix = [xyallfix; xynowfix];
    end
end
%% plot coordinates of the segmented cells (after fixing)
figure(2), hold on; imshow(t(:,:,1),[0 5e3]); hold on;
plot(xyall_rot_trans(:,1),xyall_rot_trans(:,2),'c.');
figure(1), hold on;
plot(xyallfix(:,1),xyallfix(:,2),'r.');

%xyall rotated and translated peaks (from actual live cell data, last time point)
%xyallfix peaks, moved into the frame of the montge image
%% group the new peaks into colonies (fixed data)
% 
paramfile = 'setUserParamLiveImagingAN_40Xfixed';
run(paramfile);

for k=77%51:77
outfile = [ num2str(positions2(k)) strdirfix ];
 disp(outfile);
 
load(outfile,'peaks','colonies');
clear colonies
nc = size(peaks{1},1);
pts = zeros(nc,2);
for ii=1:nc
    pts(ii,:) = (peaks{1}(ii,1:2));
end
allinds=NewColoniesAW(pts);
ngroups = max(allinds);
%ncells = cell(ngroups,1);
%Make colony structure for the single cell algorythm
for ii=1:ngroups;
    cellstouse=allinds==ii;
    colonies(ii)=dynColony(peaks{1}(cellstouse,:));
        
end

save(outfile,'colonies','-append');
end
%% check the colony
%
dirr = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-09-20-nucmasks_backtocoord/Segmented_backtoCoordTile1';
outfile = '76_testBacktocoord';
colormap = colorcube;
load(outfile,'peaks','colonies','imgfiles')
Lnuc = imgfiles(1).NucMask;

figure(2), hold on ,imshow(Lnuc,[])
for k=1:size(colonies,2)
figure(1), hold on, plot(colonies(k).cells(:,1),colonies(k).cells(:,2),'.','color',colormap(k+5,:,:),'markersize',20);
hold on
end
clear colonies
%% plot colonies onto montage
figure(2), hold on
colormap = prism;
str = '_testBacktocoord';
for k=1:77
outfile = [ num2str(k-1) str ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    
figure(2), hold on, plot(colonies(j).cells(:,1)+ac2(k).absinds(2),colonies(j).cells(:,2)+ac2(k).absinds(1),'.','color',colormap(j+8,:,:));
hold on
end
end
end
% bsxfun(@plus,peaks{1}(:,1:2),ac2(k).absinds([2 1]));


%% translate and rotate the old colonies as a whole , plot immediately onto montage

% direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
% ff = dir(direc);
% positions = (0:39);
%  load('align.mat')
%  
%  load('/Volumes/data2/tform.mat');
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [934-815, 1660-1330];%empirically determined from image
%toTranslate2 = 4*toTranslate-[377 424]; %correction empirically determined. [200 300]
toTranslate2 = [-140, 40];
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);
colormap = cool;
%str = '_lastframecol';

direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
ff = dir(direc);
positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 0;
l = 100;
for k=1:size(positions,2)
  outfile = [ num2str(k-1) strdir ];
  load(outfile,'colonies')
  if ~isempty(colonies)
    for j=1:size(colonies,2)
        for h=1:size(colonies(j).cells,2)
        if colonies(j).cells(h).onframes(end) == (l-N)
    colonies1 = bsxfun(@plus,colonies(j).cells(h).position(end,:),ac(k).absinds([2 1]));
    colonies2 = bsxfun(@minus,colonies1,mean(xyall));% ?xyall g
    colonies3 = colonies2*rotmat;
    colonies4 = bsxfun(@plus,colonies3,midpoint([2 1]));
    colonies5 = bsxfun(@plus,colonies4,toTranslate2);
    colonies(j).cells(h).position(end,:) = colonies5;
    figure(3), hold on
    plot(colonies(j).cells(h).position(end,1),colonies(j).cells(h).position(end,2),'c.','markersize',14);%'color',colormap(j+10,:,:)
        end
        
        end
    end
  end
end
    
        
% 
% % bsxfun(@plus,peaks{1}(:,1:2),ac2(k).absinds([2 1]));
% 
% xyall_rot = bsxfun(@minus,xyall, mean(xyall)); % mean subtract;
% xyall_rot = xyall_rot*rotmat; %rotate about center
% xyall_rot = bsxfun(@plus,xyall_rot,midpoint([2 1])); %add back
% 
% xyall_rot_trans = bsxfun(@plus,xyall_rot,toTranslate2);
% 
% figure; imshow(t(:,:,2),[0 5e3]); hold on;
% plot(xyall_rot(:,1),xyall_rot(:,2),'r.');
% plot(xyall(:,1),xyall(:,2),'g.');
% plot(xyall_rot_trans(:,1),xyall_rot_trans(:,2),'c.');

%% plot fixed colonies onto montage
figure(1), hold on
colormap = colorcube;
str = '_testBacktocoord';
for k=1:77
outfile = [ num2str(k-1) str ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    
figure(1), hold on, plot(colonies(j).cells(:,1)+ac2(k).absinds(2),colonies(j).cells(:,2)+ac2(k).absinds(1),'.','color',colormap(j+8,:,:));
hold on
end
end
end

%% match colonies
dirfix = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-09-20-nucmasks_backtocoord/Segmented_backtoCoordTile1');
dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/align.mat');
load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/tform.mat');
colormap = prism;
t = imread('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/test4mid.tif');
figure(1), imshow(t(:,:,1),[0 1900]);
fluordata = [];
fluordata1 = [];
fluordata2 = [];
fluordata3 = [];
allcell = [];
strfix = '_testBacktocoord.mat';
colfixall = [];% will have the centers of the colonies for the fixed data
for k=1:77
outfile = [ dirfix '/' num2str(k-1) strfix ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    colfix(:,1) = mean(colonies(j).cells(:,1)+ac2(k).absinds(2)); % mean of the x coord if the clony
    colfix(:,2) = mean(colonies(j).cells(:,2)+ac2(k).absinds(1)); % mean of y
    colSZ = size(colonies(j).cells(:,1),1);
    fluordata1 = mean(colonies(j).cells(:,5));% RFP h2b
    fluordata2 = mean(colonies(j).cells(:,6));% signaling
    fluordata3 = mean(colonies(j).cells(:,7));% cdx2
    colfixall = [colfixall;colfix]; % put data in the matrix
    fluordata = [fluordata;fluordata1 fluordata2 fluordata3 colSZ];
    %allcell = [allcell;colonies(j).cells(:,1)+ac2(k).absinds(2); colonies(j).cells(:,2)+ac2(k).absinds(1)];
%figure(1), hold on, plot(colonies(j).cells(:,1)+ac2(k).absinds(2),colonies(j).cells(:,2)+ac2(k).absinds(1),'.','color',colormap(j+8,:,:));hold on
figure(1), hold on, plot(colfix(:,1),colfix(:,2),'.','color',colormap(j+8,:,:),'markersize',10);hold on

end
end
end
% get the colonies from live cell data in the form of the matrix 
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [934-815, 1660-1330];%empirically determined from image
%toTranslate2 = 4*toTranslate-[377 424]; %correction empirically determined. [200 300]
toTranslate2 = [-140, 40];
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);
colormap = prism;
positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 0;
l = 100;
colliveall = [];% will have the centers of the colonies for the live data
colshifted = [];
col1 = [];
col2 = [];
nc = [];
for k=1:size(positions,2)
    outfile = [ dirlive '/' num2str(k-1) strdir ];
    load(outfile,'colonies')
    if ~isempty(colonies)
        for j=1:size(colonies,2)
            for h=1:size(colonies(j).cells,2)
                if colonies(j).cells(h).onframes(end) == (l-N)
                    colSZ = colonies(j).numOfCells(l-N);
                    colonies1 = bsxfun(@plus,colonies(j).cells(h).position(end,:),ac(k).absinds([2 1]));
                    colonies2 = bsxfun(@minus,colonies1,mean(xyall));% ?xyall g
                    colonies3 = colonies2*rotmat;
                    colonies4 = bsxfun(@plus,colonies3,midpoint([2 1]));
                    colonies5 = bsxfun(@plus,colonies4,toTranslate2);
                    colonies(j).cells(h).position(end,:) = colonies5;
                    %disp(colonies(j).cells(h).position(end,:));
                    
                    figure(1), hold on
                    plot(colonies5(:,1),colonies5(:,2),'.','color',colormap(j+10,:,:),'markersize',10);%
                    col1 = [col1; colonies5]; % the coordinates of all cells that were identified in the last time point
                    col2 = [col2; colSZ];     % the vector with the colony size that the cell belongs to
                end
                
            end
             
        end
        
    end
    
    
end
% and the average colony coordinates
% 
toaver = cat(2,col1,col2);      
m = max(toaver(:,3));
colaverlive = zeros(length(toaver),2); % the array with the colony center coordinates 

for k=1:m
  [r,c]= find(toaver(:,3) == k);
  if k == 1
   colaverlive(r,1:2) = toaver(r,1:2); % one cell colonies centroid coordinates
  end
   if k > 1
       for jj=1:k:(length(r))
   colaverlive(r(jj),1) = mean(toaver(r(jj):r(jj)+k-1,1));
   colaverlive(r(jj),2) = mean(toaver(r(jj):r(jj)+k-1,2)); 
      
       end
   end
   
end  
%% show the image and plot colony centroids on it from both fixed and live datasets
s = 2;% figure number
figure(s), imshow(t(:,:,2),[0 1900]); % t(:,:,2), the fixed data
load('colCenterfixed.mat')
for k=1:size(colfixall,1)
figure(s), hold on,plot(colfixall(k,1),colfixall(k,2),'.r','markersize',14);
end
load('colCenterlive.mat')
for k=1:size(colaverlive,1)
if ~isempty(colaverlive(k,:))
figure(s), hold on
plot(colaverlive(k,1),colaverlive(k,2),'.c','markersize',14);
end
end

%% match colonies based on the distance between their centroids in live vs fixed dataset

load('colCenterfixed.mat','colfixall','fluordata');
load('colCenterlive.mat','colaverlive');

% [r,c] = find(colaverlive == 0);
% colaverlive(r,:)=[];
%colfixall vs colaverlive

dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');

% load the live cell mat file 
% get the coordinates of the sifted colony
% find the match of this colony centroid to the one in the fixed data
% returm the cdx2 levels in that colony

positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 0;
l = 100;
col1 = [];
col2 = [];
tomatch = [];
for k=1:size(positions,2)
    outfile = [ dirlive '/' num2str(k-1) strdir ];
    load(outfile,'colonies2')
    if ~isempty(colonies2)
        for j=1:size(colonies2,2)
            for h=1:size(colonies2(j).cells,2)
                if colonies2(j).cells(h).onframes(end) == (l-N)
                    colSZ = colonies2(j).numOfCells(l-N);
                    
                    colonies5(h,1) = colonies2(j).cells(h).position(end,1);
                    colonies5(h,2) = colonies2(j).cells(h).position(end,2);
                    col1 = colonies5; % the coordinates of all cells that were identified in the last time point
                    col2 = colSZ;     % the vector with the colony size that the cell belongs to
                end
                
            end
     
            tomatch(1,1) = mean(nonzeros(col1(:,1))); % get the colony centroid, to test the match
            tomatch(1,2) = mean(nonzeros(col1(:,2)));
            
            if find(colaverlive(1,:)== tomatch(1,1)) == 1 % make sure that this colony centroid was identified before
                
                % here find the matching centroid in the fixed data
                d = ipdm(tomatch,colfixall,'Subset','NearestNeighbor','Result','Structure');%disp(d);
                meancdx2incol = fluordata(d.columnindex,3)/fluordata(d.columnindex,1);
                meansignalincol = fluordata(d.columnindex,2)/fluordata(d.columnindex,1);
            
            % here need to save this back into the colonies2, so that can display onthe
            % trace (at last time point)
            for h=1:size(colonies2(j).cells,2)
            colonies2(j).cells(h).fluorData(:,end) = meancdx2incol;
            end
            end
            
        end%
        
    end
    %save(outfile,'colonies2','-append');
end
