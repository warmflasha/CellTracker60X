
% to get the colonies and put them in plate1 object
function [colonies] = GetMicroColonies(matfile)

pp = load(matfile,'peaks','NucMasks');
peaks=pp.peaks;
dims = size(peaks);


for j=4:dims(2)
    
colonies=peaksToMicroColoniesAN(matfile,j);
%plate1=plate(colonies,dims,[],ff.w,[],[], matfile);
end



end