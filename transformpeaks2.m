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

strfix = '_testBacktocoord';
colfix = [];% will have the centers of the colonies for the fixed data
for k=1:77
outfile = [ num2str(k-1) strfix ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    colfix(:,1) = mean(colonies(j).cells(:,1)+ac2(k).absinds(2)); % mean of the x coord if the clony
    colfix(:,2) = mean(colonies(j).cells(:,2)+ac2(k).absinds(1)); % mean of y
    
figure(1), hold on, plot(colonies(j).cells(:,1)+ac2(k).absinds(2),colonies(j).cells(:,2)+ac2(k).absinds(1),'.','color',colormap(j+8,:,:));hold on
figure(1), hold on, plot(colfix(:,1),colfix(:,2),'*','color',colormap(j+8,:,:));hold on

end
end
end
% gert the centroids of the live colonies
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [934-815, 1660-1330];%empirically determined from image
%toTranslate2 = 4*toTranslate-[377 424]; %correction empirically determined. [200 300]
toTranslate2 = [-140, 40];
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);
colormap = cool;
positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 0;
l = 100;
collive = [];%
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
    disp(j);
    colonies(j).cells(h).position(end,:) = colonies5;
    disp(colonies(j).cells(h).position(end,:));
    
    figure(1), hold on
    plot(colonies(j).cells(h).position(end,1),colonies(j).cells(h).position(end,2),'c.','markersize',14);%'color',colormap(j+10,:,:)
        end
        
        end
         %plot(mean(colonies(j).cells(h).position(:,1)),colonies(j).cells(h).position(end,2),'c.','markersize',14);
    end
  end
end

