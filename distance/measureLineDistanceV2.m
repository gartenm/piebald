function measureLineDistanceV2(Ch1,Ch2,scale,imageNumber,fn)

%addpath 'H:\OneDrive - National Institutes of Health\projects\Piebald
%PVM\mapping' to load the OpenTifFile function

% if ~exist('fn','var')
%     fn = 'NH38-12.1.90_017.tif';
% end
%scale = 0.746540; %nm/px
% data = bfopen(fn);
% 
% ht = data{1,2};
% %scale = ht.get('Global Spacing'); %make sure to set the scale in imageJ
% scale = 1/ht.get('Global XResolution');
% 
% 
% lineCh1 = 2;
% lineCh2 = 3;
% 
% Ch1 = data{1,1}{lineCh1,1} > 1;
% Ch2 = data{1,1}{lineCh2,1} > 1;


figure(436)

imagesc(Ch1+(Ch2.*2));

NCh1 = sum(sum(Ch1))
NCh2 = sum(sum(Ch2));

% Ch1Coord = NaN(NCh1,2);
% Ch2Coord = NaN(NCh2,2);
% Ch1count = 0;
% Ch2count = 0;
% 
% for i = 1:size(Ch1,1)
%     
%     for k = 1:size(Ch1,2)
%         if Ch1(i,k)>0
%             Ch1count = Ch1count+1;
%             Ch1Coord(Ch1count,:) = [i,k];
%         end
%         if Ch2(i,k)>0
%             Ch2count = Ch2count+1;
%             Ch2Coord(Ch2count,:) = [i,k];
%         end
%     end
% end

[row,col] = find(Ch1);
Ch1Coord = [row,col];
[row,col] = find(Ch2);
Ch2Coord = [row,col];


DistMat = NaN(NCh2,1);
Ch1Ch2Dist = NaN(NCh1,6);
for i = 1:NCh1
    i
    ya = Ch1Coord(i,1);
    xa = Ch1Coord(i,2);
    for k = 1:NCh2
        yb = Ch2Coord(k,1);
        xb = Ch2Coord(k,2);
        DistMat(k) = sqrt((xb-xa)^2+(yb-ya)^2);
    end
    [dist,idx] = min(DistMat);
    Ch1Ch2Dist(i,:) = [xa,ya,Ch2Coord(idx,2),Ch2Coord(idx,1),dist,dist*scale];
end

figure(234);
hist(Ch1Ch2Dist(:,6),30);

savename = strcat(fn(1:end-3),num2str(imageNumber),'_PVMPPMdist.mat');
save(savename,'Ch1Ch2Dist','scale');
        
        


