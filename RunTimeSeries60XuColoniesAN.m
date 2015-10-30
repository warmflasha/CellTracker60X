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
for k=1:time
    
[outdatnuc,outdatcyto,Lnuc,Lcytofin] = IlastikplusWatershed_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,flag);

% here need to process the outdatnuc and outdatcyto before putting into
% peaks
% possibly empty the peaks array if it's length (number of found objects) does not match the length of nuc peaks 
peaksnuc{k} = outdatnuc;
peakscyto{k} = outdatcyto;
if size(outdatnuc,1) == size(outdatcyto,1)
peaks{k} = [outdatnuc outdatcyto];
end
NucMasks{k} = Lnuc;
CytoMasks{k} = Lcytofin;

end
dims = size(peaks);

 %save('OutfileFrame1','peaks','NucMasks','CytoMasks','dims');

end

