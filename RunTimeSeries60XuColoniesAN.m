function [peaks,dims,NucMasks,CytoMasks] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

info = h5info(ilastikfile);
info.Datasets;   

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,3);   

peaks = cell(1,time); 
peakscyto = cell(1,time);
NucMasks = cell(1,time);
CytoMasks = cell(1,time);
for k=3:time
    
[datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,flag);

peaks{k} = datacell;

%colonies=peaksToMicroColoniesAN(matfile,k);
%plate1=plate(colonies,dims,[],ff.w,[],[], matfile);

NucMasks{k} = Lnuc;
CytoMasks{k} = Lcytofin;

end
dims = size(peaks);

 %save('OutfileFrame10','peaks','NucMasks','CytoMasks','dims');

end

