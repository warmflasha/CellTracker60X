function ilastikuserparam

global userparam
%
userparam.nzslices = 21;
%%
% parameters for function primary filter
%  % number of zslices in the image
userparam.logfilter = 10;
userparam.bthreshfilter = 0.25;
userparam.diskfilter = 3;
userparam.area1filter = 100;

%%
% parameters for function secondary filter
userparam.minstartobj = 4;
userparam.minsolidity = [0.9, 0.8];
userparam.area2filter = 300;

%%
% tracing objects
% zmatch: zslices across which an object needs to be tracked. 
% if you think that one object can be present in a maximum of 5 zslices, set that as the limit.
% matchdistance: minimum distance in pixels between the centroid of two
% objectsin two different zslices to be considered as the same object.

userparam.zmatch = 21; 
userparam.matchdistance = 15;

%%
% overlap filter
userparam.overlapthresh = 60; 
userparam.imviews = 0;
%%
%  
userparam.channels = [0]; %fluorescent channel for which mRNA's are counted and need to be assigned
userparam.cmcenter = 60;
userparam.negativecontrol = 1;

end