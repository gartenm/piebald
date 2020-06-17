function batchMap12

fnlist = dir('*.czi');



N = length(fnlist);
% parfor i = 1:N
%     i
%     fnlist(i).name
% end
%     fnlist(i).name
%tic
for i = 1:N
    
    name = fnlist(i).name;
    %tic
    
    namebase = name(1:end-4);
    chkname = strcat(namebase,'__map_mollweide.tif');
    if ~isfile(chkname)
        i
        name
        piebaldmapcolocMultiChV12_12(name);

        %t1 = toc
        %tic
        savename3d = strcat(name(1:length(name)-4),'_1_sphere.mat');
        MergeProject3D(savename3d);
        %t2 = toc
        %savename3d = strcat(fn(1:length(name)-4),'1_sphere.mat');
        %MergeProject3D(savename3d);
    end
end
%totalT = toc
% for i = 1:N
%     name = fnlist(i).name;
%     savename3d(i).name = strcat(fn(1:length(name)-4),'1_sphere.mat');
% end

% parfor i = 1:N
%     i
%     %name = savename3d(i).name
%     name = fnlist(i).name
%     %piebaldmapcolocMultiChV9(name);
% 
% end