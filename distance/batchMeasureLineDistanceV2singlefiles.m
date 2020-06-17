function batchMeasureLineDistanceV2singlefiles

fnlist = dir('*.tif');

for i = 1:length(fnlist)
    batchMeasureLineDistanceStackTiff(fnlist(i).name);
    %measureLineDistanceV2(fnlist(i).name);
end