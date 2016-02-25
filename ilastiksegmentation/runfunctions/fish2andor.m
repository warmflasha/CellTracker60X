%%
function fish2andor(inputfiledir, outputdir, sample)

% from fish to andor format for further image alignment
%inputfiledir = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/images00';
%outputdir = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/mrnapercells';


for sample = 1:sample
    
    prefix = strcat(sprintf('fish%d', sample));
   
    files2move = strcat(inputfiledir, '/', strcat(prefix, '_*.tif'));
    result2move = strcat(outputdir, '/', sprintf('sample%02dout.mat', sample));
    
    newdirec = strcat(inputfiledir ,'/', prefix);
    mkdir (newdirec);
    
    movefile(files2move, newdirec);
    movefile(result2move, newdirec);
        
end

