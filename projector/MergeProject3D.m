function MergeProject3D(fn1)

if ~exist('fn1','var')
    fn1 = 'throph Coloc AS_Airyscan Processing_1_sphere.mat';
end

% 
% fnlist = dir('*1_sphere.mat');
basename = fn1(1:(length(fn1)-12));
fn2 = strcat(basename,'2_sphere.mat');

%[mapNorth,mapSouth, mapTif] = plotPol3d(fn);

savename = strcat(basename,'_map_mollweide.tif');
savename2 = strcat(basename,'_map_pol1map.tif');
savename3 = strcat(basename,'_map_pol2map.tif');


[mapNorth,mapSouth, mapTif] = plotPol3d(fn1);
imwrite(uint16(mapTif), savename, 'Compression','none');
imwrite(uint16(mapNorth), savename2, 'Compression','none');
imwrite(uint16(mapSouth), savename3, 'Compression','none');

if isfile(fn2)
[mapNorth,mapSouth, mapTif] = plotPol3d(fn2);
imwrite(uint16(mapTif), savename, 'WriteMode', 'append',  'Compression','none');
imwrite(uint16(mapNorth), savename2, 'WriteMode', 'append',  'Compression','none');
imwrite(uint16(mapSouth), savename3, 'WriteMode', 'append',  'Compression','none');
end