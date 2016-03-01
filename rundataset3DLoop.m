
function rundataset3DLoop(ilastikdirnuc,ilastikdircyto,imagedir,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto)

% Loop over seprate positions (spatial, not time points)
% loop over time points is within rundataset3D
% see rundataset3D for input parameters

ff=readAndorDirectory(imagedir);
pp = max(ff.p);
positions = 0:(positions-1);

for jj=1:pp
pos = position(jj);

outfile = ([ num2str(pos) '_' num2str(outfile)]);


rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto)
end
