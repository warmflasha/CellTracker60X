
% to get the colonies and put them in plate1 object
function [plate1] = GetMicroColonies(matfile)

pp = load(matfile,'peaks','NucMasks','dims');
peaks=pp.peaks;
dims = size(peaks);


for j=1:dims(2)
    
colonies(j)=peaksToMicroColoniesAN(matfile,j);
end
plate1=plate(colonies(j),dims,[],ff.w,[],[], matfile);


end