function nucleimrnacheck(masterCC, inuc, nucleilist, stats, objno, cellmrna, zrange)

zstart = zrange(1);
zend = zrange(end);

nucleiinfo = cell(1,size(nucleilist,1));

for i = 1:size(nucleiinfo,2)
    
    object = nucleilist(i,:);
    objectmatchcol = find(~isnan(object));
    
    objectmatchz = objectmatchcol + zstart -1;
    objectsmatched = nucleilist(i,objectmatchcol);
    
    for match = 1:numel(objectsmatched)
        centroidmatch = stats{objectmatchz(match)}(objectsmatched(match)).Centroid;
        nucleiinfo{i}(match,:) = [centroidmatch objectmatchz(match)];
    end
    
    xcenter = 0.5*(max(nucleiinfo{i}(:,1))+min(nucleiinfo{i}(:,1)));
    ycenter = 0.5*(max(nucleiinfo{i}(:,2))+min(nucleiinfo{i}(:,2)));
    zcenter = 0.5*(max(nucleiinfo{i}(:,3))+min(nucleiinfo{i}(:,3)));
    
    nucleicenter(i,:) = [xcenter, ycenter, zcenter];
end

newcenterlist = nucleicenter;


masterLntr = labelmatrix(masterCC);
rgbLntr = label2rgb(masterLntr, 'jet', 'k', 'shuffle');
overlayntr = zeros(size(rgbLntr));

zstart = zrange(1);
zend = zrange(end);
for c1 = 1:3
    overlayntr(:,:,c1) = 0.5*mat2gray(rgbLntr(:,:,c1)) +...
        mat2gray(sum(inuc(:,:,zstart:zend),3));
end
figure; imshow(overlayntr);

hold on;

for i = 1:size(newcenterlist,1)
    text(newcenterlist(i,1), newcenterlist(i,2), int2str(i), 'Color', 'w', 'FontSize', 12);
    
end

figure; imshow(overlayntr);
hold on;

for i = 1:size(newcenterlist,1)
    text(newcenterlist(i,1), newcenterlist(i,2), int2str(i), 'Color', 'w', 'FontSize', 12);
    
end



for i = 1:numel(objno)
    hold on;
    plot(cellmrna{objno(i)}{1}(:,3), cellmrna{objno(i)}{1}(:,4), 'g.');
end


