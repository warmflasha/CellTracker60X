function [smoothtrace] = smoothVAR(vect,N)

smoothwindow = N;
dat = power(vect,2);
a = find(isnan(dat));
dat(a) = 0;
smoothtrace = imfilter(dat,ones(smoothwindow,1))/smoothwindow;

end