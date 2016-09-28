function [outdat,Lnuc] = runmaskoneANdataonlynucsegm(pnuc,inuc,paramfile,paramfile3D,imagedir,pl,pos,timegroup)


eval(paramfile);
eval(paramfile3D);
global userParam;

% sapna code tracking 
pmasks = primaryfilterAN(pnuc,userParam.probthresh_nuc, userParam.area1filter);

% zrange: where the nuclei are in z
[zrange] = secondaryfilterAN(pmasks, userParam.minstartobj);
smasks = pmasks;
%[zrange, smasks] = secondaryfilterAN(pmasks, userParam.minstartobj, userParam.minsolidity);
 if zrange == 0;
    outdat = [];
    Lnuc = pmasks(:,:,1);
    disp('zrange is zero1');
    return
 end

if userParam.flag ==1
    for k=1:size(smasks,3) 
        [~,smasks(:,:,k)] = UnmergetwonucleiGeneral(smasks(:,:,k));
      
    end
        
end

if userParam.flag ==0
    disp('no unmerge')
end
% get the acual tracking nucleilist = tracked objects, labled 1-N,CC -
% pixelidxlist of all objects tracked in all planes
[PILsn, PILsSourcen, CC, masterCCn, stats, nuclein1, zrange] = traceobjectszdistinct(smasks, userParam.matchdistance, zrange, size(zrange,2));%size(zrange,2)userParam.zmatch

if ~iscell(CC) 
  pmasks = primaryfilterAN(pnuc,userParam.probthresh_nuc, userParam.area1filter);
  [zrange, smasks] = secondaryfilterAN(pmasks, userParam.minstartobj);
  if zrange == 0;
    outdat = [];
    Lnuc = pmasks(:,:,1);
   
    disp('zrange is zero2');
    return
end
  [PILsn, PILsSourcen, CC, masterCCn, stats, nuclein1, zrange] = traceobjectszdistinct(smasks, userParam.matchdistance, zrange, size(zrange,2));%size(zrange,2)userParam.zmatch  
     
end
if zrange == 0;
    outdat = [];
    Lnuc = smasks(:,:,1);
    
    disp('zrange is zero3');
    return
end
 % use nucleilist to relabel the tracked objects with unique labels
[newmask_lbl] = lblmask_3Dnuc(CC,nuclein1);

% here can select the good plane, make different var pl and set it to the numeric value of the zslice that is the best   

if size(pnuc,3)>1
    selectZ = [2 3];
else
    selectZ = 3;
end
 [datacell,Lnuc] = nucOnlyIlastik2peaks_3Dsegm(newmask_lbl,inuc,paramfile,imagedir,selectZ,pos,timegroup);
 outdat = datacell;
 
 
end
