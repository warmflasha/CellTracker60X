%function montagetomultxy(filep)

filep = 'an4.xyz';
%%
%transforming montage to multi-XY scan text file.
% reading file

file2 = fopen(filep);
file2r = textscan(file2, '%s');

% finding positions of relevant text

n1 = size(file2r{1},1);
posnote = [{'XYFields'}; {'MontageOffsets'}; {'MF'}; {'AutofocusPositions'}];
%posnote = [{'XYFields'}; {'MontageOffsets'}; {'MF'}];

lastpart = [{'ZScan='}; {'Montage='}; {'Plates='}; {'MontageOffsets='}; {'MF Scan Present=1'}];
%lastpart = [{'AutofocusPositions='}; {'ZScan='}; {'Montage='}; {'Plates='}; {'MontageOffsets='}; {'MF Scan Present=1'} ];


for i = 1:size(posnote,1)
    for j = 1:n1
      txt = file2r{1}{j};
      if(~isempty(strfind(txt,posnote{i})))
          pos(i) = j;
          break;
      end
    end        
end
%%
%transfer all montage fields to a cell array

ffield = pos(2) + 1;
lfield = pos(3) - 1;
m  =1;

for i = ffield:lfield
    newfields{m} = file2r{1}{i};
    m = m+1;
end

%%
% adding initial value to newfileds
sfield = pos(1) + 1; 
startval = file2r{1}{sfield};
startval = str2num(startval);

%%
for i = 1:size(newfields,2)
    val = str2num(newfields{i});
    valn = val + startval;
    for j = 1:numel(valn)
        a{j} = num2str(valn(j));
    end
      
     vals = strcat(a{1}, ',', a{2}, ',', a{3});
     newfields{i} = vals;
    
end

%%
nfields = m-1;

%%
% modify the no. of fields to nfields

str1 = file2r{1}{pos(1)};
i1 = strfind(str1, '=');
str2 = num2str(nfields);

str3 = strcat(str1(1:i1), str2);

%new autofocus option
str11 = file2r{1}{pos(4)};
i11 = strfind(str11, '=');
str4 = strcat(str11(1:i11), str2);

%%
newfilename = 'multiplexy_AN.xyz';

% writing a new file

    n1 = size(file2r{1},1);
    newfileid = fopen(newfilename, 'w');
    
   break1 = pos(1)-1;
   
   for i = 1:break1  % initial part
     fprintf(newfileid, '%s\n', file2r{1}{i,:});
   end
     
     fprintf(newfileid, '%s\t', str3); % new XY fields
     
   for m = 1:nfields-1
     fprintf(newfileid, '%s\t', newfields{m});
   end
   
     fprintf(newfileid, '%s', newfields{nfields});
     
     fprintf(newfileid, '\n%s ', str4);  % new Autofocus
  %%   
     
  af = num2str(8);
     for m = 1:nfields-1
         fprintf(newfileid, '%s\t', af);
     end
%    
     
         fprintf(newfileid, '%s', af);
     

     for i = 1:size(lastpart,1)
         fprintf(newfileid, '\n%s', lastpart{i});
     end
     
      fprintf(newfileid, '\n' );
    
    %%
    fclose(newfileid);
    
    
        
    
