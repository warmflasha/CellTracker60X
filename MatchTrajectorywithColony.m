% plot the signaling value on top of the image
function MatchTrajectorywithColony(n1,n2,m,N,matfile,coltoplot)
% frb2016 data
ilastikdircyto = ('/Users/warmflashlab/Desktop/Feb2016ilastik_CytoMasks');
ilastikdirnuc = ('/Users/warmflashlab/Desktop/Feb2016ilastik_NucMasks');
imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Nuc_raw_data');%rawimages_nuc
imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Cyto_raw_data');%rawimages_cyto

% tiling data July 7
%  ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset');%/cyto_mask_z1.h5'
% ilastikdirnuc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset');%/nuc_mask_z1.h5'
% imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/nuc_raw');%rawimages_nuc%
% imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/cyto_raw');%rawimages_cyto


%imagedir =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/testPos1_rawImgs');
%pos = 12;
%tpt =1;

%timegroup = 1;
timegroup = [];

chan = [1 2];
chanal = 1;

  positions = [0 2 3 4 5 8 9 11 12 13 14 15 16 17 18 20 21 22 25 26 32 33 10]; %  february 3 dataset
 % positions = [15]; %  february 3 dataset
pl = 1;%3
% strcyto = 'cytomask3DFebset'; %CytoMasks3Djan8set  %cytomask3DFebset
% strnuc = 'nucmask3DFebset';   %NucMasks3Djan8set   %nucmask3DFebset
strcyto = 'cyto_mask_'; 
strnuc = 'nuc_mask_';
pos = positions(N);
[ilastikCytoAll] = FindPositionMasks(ilastikdircyto,pl,pos,strcyto);    % get the specific position ilastik masks (all z projections)
[ilastikNucAll] = FindPositionMasks(ilastikdirnuc,pl,pos,strnuc);

imgfilescyto = struct;
imgfiles = struct;

% read raw images
load(matfile);
[imgsnuc_reader]   =  getrawimgfiles(imagedir1,pl, pos,timegroup,chanal(1));        % get the raw images for that position and merge them into a 3d format
[imgscyto_reader] =   getrawimgfiles(imagedir2,pl, pos,timegroup,chanal(1));
nT = imgsnuc_reader{1}.getSizeT;  

numcol = size(colonies,2);% how many colonies were grouped within the image
% colors = colormap(colorcube);
% cmap = colors(1:20,:);close all
c = {'r','g','b','r','c','m','g','r','g','b','r','c','m','g'};

if coltoplot > numcol
    disp('there are no more colonies in this image')
    return
end
%m = 2;% which slice to take
for k=n1:n2
    
        curr = imgsnuc_reader{m}.getIndex(0,0, k - 1) + 1;
        img = bfGetPlane(imgsnuc_reader{m},curr);
        
        %for j=coltoplot
        for i=1:size(colonies(coltoplot).cells,2)
        figure(k),subplot(1,2,1),imshow(img,[]);hold on
        subplot(1,2,1),plot(colonies(coltoplot).cells(i).position(:,1),colonies(coltoplot).cells(i).position(:,2),'*','color',c{i});hold on%cmap(i,:)
        %hold on, text(colonies(j).cells(i).fluorData(k,2)./colonies(j).cells(i).fluorData(k,3)+10,num2str(colonies(j).cells(i).fluorData(k,2)./colonies(j).cells(i).fluorData(k,3)),'Color','m');
        subplot(1,2,2),plot(colonies(coltoplot).cells(i).onframes',colonies(coltoplot).cells(i).fluorData(:,2)./colonies(coltoplot).cells(i).fluorData(:,3),'-*','color',c{i});hold on

    end
    %end
    
end
ylim([0 2.5]);
end