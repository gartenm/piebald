function batchMeasureLineDistanceStackTiff(fn)

%fnlist = dir('*.tif');
if ~exist('fn','var')
    fn = '20191220 GDB132 A4 PfNCR1 control.tif';
end
data = bfopen(fn);

ht = data{1,2};
%scale = ht.get('Global Spacing'); %make sure to set the scale in imageJ
scale = 1/ht.get('Global XResolution');

N = size(data{1,1},1)/3;

for i=1:N
    idx = (i-1)*3;
    
    lineCh1 = idx+2;
    lineCh2 = idx+3;
    Ch1 = data{1,1}{lineCh1,1} > 1;
    Ch2 = data{1,1}{lineCh2,1} > 1;

    measureLineDistanceV2(Ch1,Ch2,scale,i,fn)
    
end

