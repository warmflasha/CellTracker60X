
% test the watershed + ilastik analysis of 60X data


direccurr = ('.');
cd(direccurr);
direc = ('Z:\Anastasia\LiveCellImagingGFPs4RFPh2b\SingleCellSignalingAN_20150805_123245 PM');
%direc = ('/Volumes/data/Anastasia/LiveCellImagingGFPs4RFPh2b/SingleCellSignalingAN_20150805_123245 PM');
flag = 0;
ilastikfile = 'NucMask0001batch.h5';
ilastikfilecyto = 'CytoMask0001batch.h5';
pos = 1;
zplane = 4;
img = 1;

%to run img = 1,first image


 [datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,img,flag);


% to run all time points
% 
[peaks,dims,NucMasks,CytoMasks] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag);

