function [img,meta] = OpenCziZStack(fn,WorkChannel)


%files to load
if ~exist('fn','var')
    fn = 'cell1-Airyscan Processing-01.czi';
end

if ~exist('WorkChannel','var')
    WorkChannel = 1;
end


CZI = bfopen(fn);

metastring = CZI{1,2};

% % create example Java hashtable and add some key/value pairs
% h = java.util.Hashtable;
% h.put('Name', 'John');
% h.put('Country', 'UK');
% % retrieve the value for a specific key
% h.get('Name');
% % retrieve all key names
% allKeys = arrayfun(@char, h.keySet.toArray, 'UniformOutput', false);
% % retrieve all key values
% allValues = cellfun(@(x) h.get(x), allKeys, 'UniformOutput', false);


out = metastring.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1');
scale(1) = str2num(out);
scale(2) = str2num(out);
out = metastring.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1');
scale(3) = str2num(out);
out = metastring.get('Global Information|Image|SizeZ #1');
NZ = str2num(out);
out = metastring.get('Global Information|Image|SizeC #1');
NC = str2num(out);


% metaCommaIdx = strfind(metastring,',');
% ScalePos1 = strfind(metastring,'ScalingZ #1');
% ScalePos2 = metaCommaIdx(find(metaCommaIdx>ScalePos1,1,'first'));
% ScaleString = metastring(ScalePos1:ScalePos2);
 
%Im = CZI{1,1};
%sigma = 2; %gaussian blur


%get information on file, find number of channels
% header = char(CZI{1,1}(1,2));
% if ~isempty(strfind(header,'Z=1/'))
%     Z(1) = strfind(header,'Z=1/')+4;
%     Z(2) = strfind(header,'; C');
%     NZ = header(Z(1):Z(2))
% end
% Z(1) = strfind(header,'C=1/')+4;
%Z(2) = strfind(header,'; C')



%copy cell to matrix, convert to 32bit, double
% NCh = str2num(header(Z(1):end));
% NStack = size(CZI{1,1},1);

img = NaN(NZ,size(CZI{1,1}{1,1},1),size(CZI{1,1}{1,1},2));

for i = 1:NZ
    idx = NC*(i-1)+WorkChannel;
    img(i,:,:) = CZI{1,1}{idx,1};
end

% figure(2);
% imagesc(squeeze(img(23,:,:)))


meta.scale = scale;
