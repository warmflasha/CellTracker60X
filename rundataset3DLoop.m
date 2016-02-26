
function rundataset3DLoop(ilastikdirnuc,ilastikdircyto,imagedir,paramfile,timegroup,outfile,paramfile3D)

% this is the loop over seprate positions 
ff=readAndorDirectory(imagedir);
pp = max(ff.p);
positions = 0:(positions-1);

for jj=1:pp
pos = position(jj);

outfile = ([ num2str(pos) '_' num2str(outfile)]);

rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir,pos,paramfile,timegroup,outfile,paramfile3D)
end
