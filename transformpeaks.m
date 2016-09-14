% code to transform the paeks nto the coordinates of the montage nd then 
% find the colonies
%%
% load outfile
% load the last peaks and make an image with peaks(:,1:2) as the centroids of objects/ or even colonies; make a directory of such images named as Andor
% then make a montage out of those peaks and transform it as above)
direc2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/lastTlivecellset/');

[ac, fi]=alignManyPanelsAndorZstackMontage(direc2,[10 4],0,1);% 
ff = dir('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
[ac2, fi2]=alignManyPanelsAndorZstackMontage('good_Backtocoord40xBF2_20160709_53148 PM',[7 11],0,1);


%direc3 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/PeaksIntoMontage_lastT/');
%imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/nuc_raw');%rawimages_nuc%
t = imread('test4mid.tif');
imshow(t(:,:,2),[]); % 2 - green chanel, 1 - red(4x10)
imcontrast
clear peaksnew
pl = 2;
positions = (0:39);
strdir = '_40X_imprBGandSegm.mat';
N = 2;
%ff2 = readAndorDirectory(imagedir1);
chan = [1];
chanal = 1;
timegroup = [];
l = size(peaks,2);
    for k=1:size(positions,2) % same as the loop over image numbers
        outfile = [ num2str(positions(k)) strdir ];
        load(outfile,'peaks','imgfiles','imgfilescyto');

        peaksnew(1:size(peaks{l-N},1),1) = peaks{l-N}(:,1) + ac(k).absinds(1)*ones(size(peaks{l-N},1),1)+toTranslate2(1)*ones(size(peaks{l-N},1),1); %+toTranslate2(2) ac.absind(2) the x axis
        peaksnew(1:size(peaks{l-N},1),2) = peaks{l-N}(:,2) + ac(k).absinds(2)*ones(size(peaks{l-N},1),1)+toTranslate2(2)*ones(size(peaks{l-N},1),1);
        figure(1), hold on
        plot(peaksnew(:,1),peaksnew(:,2),'*b','markersize',12);
    end



  %%


[ac, fi]=alignManyPanelsAndorZstackMontage(direc2,[10 4],0,1);% here need to get the fluor, corresponding to the last time points 
[ac2, fi2]=alignManyPanelsAndorZstackMontage('good_BacktoCoord2_Smad4H2B_CDX2_20160709_60748 PM',[7 11],0,1);

load('tform.mat');

mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [934-815, 1660-1330];%empirically determined from image


f2 = fi2;
f1reg = imwarp(fi,mytform);
toTranslate2 = 4*toTranslate-[377 424]; %correction empirically determined. [200 300]

f1reg_trans = zeros(size(f2));
f1reg_trans(toTranslate2(2):(toTranslate2(2)+size(f1reg,1)-1),toTranslate2(1):(toTranslate2(1)+size(f1reg,2)-1)) = f1reg;

f2(size(f2,1):size(f1reg_trans,1),:) = 0;%11265:11842

img2output = cat(3,uint16(f1reg_trans),uint16(f2),zeros(size(f2)));
imwrite(img2output,'test4mid.tif');
    
    
    