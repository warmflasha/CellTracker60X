%%intensity based registration
% approach to match the trajectory with fate:
% as the after mage need to use the segmented binary mask of the stained
% chip (taken after live cell imaginf was done)

% as the image before, need to use the last timepoint of the live cell data ( the binary mask of the best zplane ) 
%
% then just register those images based on the location of the cell object
% in the 'before image'
% image registration test, tiling 2, 40 X
aft = imread('BBieldBacktoCdx2pos_tiling2_f0000toregister.tif');
bf = imread('Tiling2Coord_afterStopLiveCell_f0000.tif');

aft1 = aft<650;
aft2 = bwareafilt(aft1,[100 10000]);
aft3 = imfill(aft2,'holes');
aft4 = imdilate(aft3,strel('disk',2));

bf1 = bf<3000;
bf2 = bwareafilt(bf1,[700 10000]);
bf3 = imfill(bf2,'holes');
bf4 = imdilate(bf3,strel('disk',2));


figure, imshow(bf4,[]);
figure, imshow(aft4,[]);

figure, imshowpair(bf,aft);% 'mpntage'
%----------------object-based registration
cpselect(aft4,bf4);

trform = fitgeotrans(movingPoints, fixedPoints,'similarity');
regim = imwarp(aft, trform);
regim2 = imwarp(bf, trform);
figure, imshow(regim,[]);
%-------------------------------intensity-based registration
[optimizer,metric] = imregconfig('monomodal');
% synax for registration function:   moving_reg = imregister(moving,fixed,transformType,optimizer,metric)
moving_reg = imregister(bf,aft3,'translation',optimizer,metric);
figure, imshow(moving_reg,[])
%% specific control object in both images-based registration

outfile = '37_test2.mat';
load(outfile);
j = 83; % take the last usable time point ( to match the cells that are fixed on the chip right after live cell was stopped)
Lnuc = imgfiles(j).NucMask;
Lnuc1 = im2bw(Lnuc);

bf = Lnuc1;%imread('bf_tp82_pos3_toregister.tif');
aft = imread('aft_backtoCdx2pos_tiling2_f0036_toregister.tif');

aft1 = aft>800;
aft2 = bwareafilt(aft1,[700 7000]);
aft3 = imfill(aft2,'holes');
aft4 = imerode(aft3,strel('disk',1));

figure, imshow(bf,[]);
figure, imshow(aft3,[]);

figure, imshowpair(bf,aft3);% 'mpntage'
cpselect(aft3,bf);

trform = fitgeotrans(movingPoints1, fixedPoints1,'similarity');
regim = imwarp(aft, trform);
regim2 = imwarp(bf, trform);
figure, imshow(regim,[]);


%% intensity based registraion
% image registration test, tiling 1, 40X
bf = imread('BFiels_tileCoord_oneZ_f0000.tif');
aft = imread('backtocoord40xBF2_f0016.tif');

figure, imshow(bf,[]);
figure, imshow(aft,[]);

figure, imshowpair(bf,aft);% 'montage'

[optimizer,metric] = imregconfig('monomodal');

% synax for registration function:   moving_reg = imregister(moving,fixed,transformType,optimizer,metric)


moving_reg = imregister(bf,aft,'rigid',optimizer,metric);
figure, imshow(moving_reg,[])
%% specific control object in both images-based registration

bf = imread('BFiels_tileCoord_oneZ_f0000.tif');
aft = imread('backtocoord40xBF2_f0016.tif');

figure, imshow(bf,[]);
figure, imshow(aft,[]);
figure, imshowpair(bf,aft);
cpselect(aft,bf);

trform = fitgeotrans(movingPoints5, fixedPoints5,'similarity');
regim = imwarp(aft, trform);
figure, imshow(regim,[]);


