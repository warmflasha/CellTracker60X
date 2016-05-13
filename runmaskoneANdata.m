function [outdat,Lnuc,Lcytofin] = runmaskoneANdata(pnuc,pcyto, inuc,icyto,timegroup,paramfile,paramfile3D)


eval(paramfile);
eval(paramfile3D);
global userParam;

% sapna code tracking 
pmasks = primaryfilter(pnuc,userParam.logfilter, userParam.bthreshfilter, userParam.diskfilter, userParam.area1filter);

% here can insert the UnmergeTwoNuclei function ( after the masks are
% already binary

% zrange: where the nuclei are in z
[zrange, smasks] = secondaryfilter(pmasks, userParam.minstartobj, userParam.minsolidity, userParam.diskfilter, userParam.area2filter);
 if zrange == 0;
    outdat = [];
    Lnuc = pmasks(:,:,1);
    Lcytofin = im2bw(pcyto(:,:,1),0.9);
    
    return
 end
% here can find the max object in the image and then based on this info set
% the parameter for the minimum area to unmerge 


if userParam.flag ==1
for k=1:size(smasks,3)
         [~,smasks(:,:,k)] = UnmergetwonucleiGeneral(smasks(:,:,k));
 end
end
% do again, to make sure that if three cells were merged, they are cut
% after the second round og Unmerge2nuclei
for k=1:size(smasks,3)
         [~,smasks(:,:,k)] = UnmergetwonucleiGeneral(smasks(:,:,k));
 end

if userParam.flag ==0
    disp('no unmerge')
end
% get the acual tracking nucleilist = tracked objects, labled 1-N,CC -
% pixelidxlist of all objects tracked in all planes
[PILsn, PILsSourcen, CC, masterCCn, stats, nuclein1, zrange] = traceobjectszdistinct(smasks, userParam.matchdistance, zrange, size(zrange,2));%size(zrange,2)userParam.zmatch

if ~iscell(CC) 
  pmasks = primaryfilter(pnuc,userParam.logfilter, userParam.bthreshfilter, userParam.diskfilter, userParam.area1filter);
  [zrange, smasks] = secondaryfilter(pmasks, userParam.minstartobj, userParam.minsolidity, userParam.diskfilter, userParam.area2filter);
  if zrange == 0;
    outdat = [];
    Lnuc = pmasks(:,:,1);
    Lcytofin = im2bw(pcyto(:,:,1),0.9);
    
    return
end
  [PILsn, PILsSourcen, CC, masterCCn, stats, nuclein1, zrange] = traceobjectszdistinct(smasks, userParam.matchdistance, zrange, size(zrange,2));%size(zrange,2)userParam.zmatch  
     
end
if zrange == 0;
    outdat = [];
    Lnuc = smasks(:,:,1);
    Lcytofin = im2bw(pcyto(:,:,1),0.9);
    
    return
end
 % use nucleilist to relabel the tracked objects with unique labels
[newmask_lbl] = lblmask_3Dnuc(CC,nuclein1);

% [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userParam.overlapthresh);

% leave only planes with cells in the cyto channel too
goodk = zeros(size(zrange,2),1);
 for k=1:size(zrange,2)
 N = size(nuclein1,1);
 a =  find(isnan(nuclein1(:,k)));
 
 if size(a,1)< N 
     goodk(k,1) = k;
 end
 goodk = nonzeros(goodk);
 end
  
 for k=1:size(goodk,1)
 %pmaskscyto1(:,:,k) = im2bw(pcyto(:,:,zrange(goodk(k))),userParam.probthresh_cyto);
 pmaskscyto(:,:,k) = imfill(pcyto(:,:,zrange(goodk(k)))> userParam.probthresh_cyto,'holes');
 %pmaskscyto(:,:,k) = bwareafilt(pmaskscyto1(:,:,k),[0 50000]);           % here filter the cyto masks by area ( remove realy huge cytoplasms )
 icyto_new(:,:,k) =icyto(:,:,zrange(goodk(k)));                         % inuc, icyto = still have 5 layers instack, need to leave only the ones
 inuc_new(:,:,k) =inuc(:,:,zrange(goodk(k)));                           % with the nuclei in them ( as was determined by zrange)
 end
 
 

 % get the labeled cyto masks ( nuclei not subtracted in the
 % GetVoronoiCells3D function since need them filled in order to remove
 % cytos without nuc)
 
 [maskzcyto_lbl] = GetVoronoiCells3D(newmask_lbl,pmaskscyto);
 
 % now maskzcyto and newmask are 3D masks that need to be applied to
 % icyto,inuc, which are all now of the same dim
 % all the masks are labeled
   
 [datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks_3Dsegm(newmask_lbl,maskzcyto_lbl,inuc_new,icyto_new,paramfile,zrange);
 outdat = datacell;
 
 
end





