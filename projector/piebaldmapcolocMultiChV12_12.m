function piebaldmapcolocMultiChV12_12(datafn)

%NChannels = 3;
%LeadChannel = 1;
%load image
P3D = 1;
if ~exist('datafn','var')
    datafn = 'cell1.tif';
end

savename = strcat(datafn(1:length(datafn)-4),'_map_comp.tif');
savename2 = strcat(datafn(1:length(datafn)-4),'_map_comp_stretch.tif');
savename3 = strcat(datafn(1:length(datafn)-4),'_map_comp_stretch_scale.tif');

out = piebladmapV12(datafn,1,0,P3D);
imwrite(uint16(out.map), savename, 'Compression','none');
imwrite(uint16(out.map2), savename2, 'Compression','none');
imwrite(uint16(out.map3), savename3, 'Compression','none');

out = piebladmapV12(datafn,2,0,P3D,out.centerX,out.centerY,out.centerZ);
imwrite(uint16(out.map), savename, 'WriteMode', 'append',  'Compression','none');
imwrite(uint16(out.map2), savename2, 'WriteMode', 'append',  'Compression','none');
imwrite(uint16(out.map3), savename3, 'WriteMode', 'append',  'Compression','none');

% out = piebladmapV2(datafn,3,1,out.centerX,out.centerY,out.centerZ);
% imwrite(uint16(out.map), savename, 'WriteMode', 'append',  'Compression','none');
% imwrite(uint16(out.map2), savename2, 'WriteMode', 'append',  'Compression','none');
% imwrite(uint16(out.map3), savename3, 'WriteMode', 'append',  'Compression','none');

% for i=1:NChannels
%     if i==1
%         out = piebladmap(datafn{1})
%     else
%         out = piebladmap(datafn{1},out.centerX,out.centerY,out.centerZ,i)
%     end
% end
