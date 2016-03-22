% info on where the images are and whidh position / timepoint to run

ilastikdircyto = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/CytoMasks3D');
ilastikdirnuc = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/NucMsks3D');
imagedir1 = ('/Users/warmflashlab/Desktop/3DanalysisRAWimg_W0');%rawimages_nuc
imagedir2 = ('/Users/warmflashlab/Desktop/3DanalysisRAWimg_W1');%rawimages_cyto

%imagedir =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/testPos1_rawImgs');
%pos = 12;
%tpt =1;
timegroup = 1;
chan = [1 2];
chanal = 1;
paramfile = 'setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN';

outfile = '3D_20hr_test_xyz.mat';


positions = [1 2 5 7 8 12 14 15 18 19 20 24 25 27 29 30]; %positions that were processed

pl = 5;
strcyto = 'CytoMasks3Djan8set';
strnuc = 'NucMasks3Djan8set';
%for k=1:length(positions)
k = 14;
    pos = positions(k);
rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto,chanal);
%end