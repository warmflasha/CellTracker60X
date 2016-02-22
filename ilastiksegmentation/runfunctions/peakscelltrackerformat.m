function peaks = peakscelltrackerformat(peaks1)

ncolumn = size(peaks1,2);
nchannels = ncolumn-3;

peaks(:,1:3) = peaks1(:,1:3);


 % arbitrary values for columns 3, 4, 5 - just to make the output
 % compatible with celltracker peaks output.
 
%peaks(:,3) = 1000;
peaks(:,4) = -1;
peaks(:,5) = 100;

channelno = 4;

for i = 6:2:6+(nchannels-1)*2
    peaks(:,i) = peaks1(:,channelno);
    peaks(:,i+1) = peaks1(:,channelno);
    
    channelno = channelno+1;
end
end

    
