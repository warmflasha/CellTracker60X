function [peaks,dims,NucMasks,CytoMasks, colonies] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

info = h5info(ilastikfile);
info.Datasets;   

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,3);   

peaks = cell(1,time); 
peakscyto = cell(1,time);
NucMasks = cell(1,time);
CytoMasks = cell(1,time);
for k=1:time
    
[datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,flag);

peaks{k} = datacell;

%colonies=peaksToMicroColoniesAN(matfile,k);
%plate1=plate(colonies,dims,[],ff.w,[],[], matfile);

NucMasks{k} = Lnuc;
CytoMasks{k} = Lcytofin;

end
dims = size(peaks);

colonies=peaksToMicroColoniesAN(peaks);% fro each time frame % here the coonies is a cell array : each cell is a colony object

 %save('Outfile','peaks','NucMasks','CytoMasks','colonies');

 
end

