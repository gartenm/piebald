function saveDataToTable

fnlist = dir('*PVMPPMdist.mat');
N = length(fnlist)
outfn = 'DistanceTable.xls';
outfn2 = 'scale.xls'
datastr = {};
for i = 1:N
    datastr{i} = load(fnlist(i).name);
end

SV = ('A':'Z');

for i = 1:N
    dataset = datastr{1, i}.Ch1Ch2Dist(:,6);
    size(dataset)
    if i <= 26
        corner = strcat(SV(i),'1');
    else
        corner = strcat('A',SV(i-26),'1');
    end
    %corner = strcat(SV(i),'1');
    
    writematrix(dataset,outfn,'Range',corner);%,'WriteMode','append')
    writematrix(datastr{1, i}.scale,outfn2,'Range',corner);
end