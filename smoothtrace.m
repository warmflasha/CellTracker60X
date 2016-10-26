function [smoothtrace] = smoothtrace(vect,N)

smoothwindow = N;
dat = vect;
a = find(isnan(dat));
dat(a) = 0;
smoothtrace = imfilter(dat,ones(smoothwindow,1))/smoothwindow;

end