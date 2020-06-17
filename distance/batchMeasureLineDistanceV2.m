function batchMeasureLineDistanceV2

fnlist = dir('*.tif');

for i = 1:length(fnlist)
    measureLineDistanceV2(fnlist(i).name);
end