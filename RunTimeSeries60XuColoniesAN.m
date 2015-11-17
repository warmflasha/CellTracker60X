function [peaks,dims,imgfilescyto,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

info = h5info(ilastikfile);
info.Datasets;   

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,3);   

peaks = cell(1,time); 
peakscyto = cell(1,time);

imgfiles = struct([]);
imgfilescyto = struct([]);
for k=1:time
    
[datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,flag);

peaks{k} = datacell;
if isempty(datacell)
continue
end

imgfiles(k).compressNucMask = compressBinaryImg(Lnuc, size(Lnuc) );
imgfilescyto(k).compressNucMask = compressBinaryImg(Lcytofin, size(Lcytofin) );

%NucMasks{k} = Lnuc;
%CytoMasks{k} = Lcytofin;

end
dims = size(peaks);
%colonies=peaksToMicroColoniesAN(peaks);% for each time frame % here the colonies is a cell array : each cell is a colony object

 %save('Outfile','peaks','NucMasks','CytoMasks','colonies');

 
end

