function [img,meta] = OpenTifFile(fn,WorkChannel)

data = bfopen(fn);

%read meta data
omechar = char(data{1,4}.dumpXML());
%get scale
ValueBrakets = strfind(omechar,'"')
dimensions = ['X','Y','Z'];
for c = 1:length(dimensions)
    searchSTR = strcat('PhysicalSize', dimensions(c), '="');
    ind1 = strfind(omechar,searchSTR)+length(searchSTR);
    if ~isempty(ind1)
    ind2 = [ValueBrakets(ValueBrakets > ind1)];
    ind2 = ind2(1)-1;
    scalestring = omechar(ind1:ind2);
    scale(c)= str2double(scalestring);
    else
        scale(c)=NaN;
    end
end
%get info about the stack
sliceinfo = data{1,1}(:,2);
%ValueBrakets = strfind(sliceinfo,'/')

searchSTR = strcat('SizeC="');
ind1 = strfind(omechar,searchSTR)+length(searchSTR);
ind2 = [ValueBrakets(ValueBrakets > ind1)];
ind2 = ind2(1)-1;
scalestring = omechar(ind1:ind2);
%scale(c)= str2double(scalestring);

NChannels = str2double(scalestring);

if NChannels > 1
    dimensions = ['Z','C'];
    %dimensions = ['Z'];

    %find images from 1st channel
    for c = 1:length(dimensions)
        for i = 1:length(sliceinfo)
            scalestring = char(sliceinfo(i));
            ValueBrakets = strfind(scalestring,'/')
            searchSTR = strcat(dimensions(c), '?=');
            
            ind1 = strfind(scalestring,searchSTR)+length(searchSTR);
            if ~isempty(ind1)
            ind2 = [ValueBrakets(ValueBrakets > ind1)];
            ind2 = ind2(1)-1;
            %scalestring = sliceinfo(i);
            scalestring = scalestring(ind1:ind2);
            imageinfo(i,c)= str2double(scalestring);
            else
                imageinfo(i,c)=NaN;
            end
        end
    end
    imageinfo
    img = data{1,1}(:,1);

    for i = 1:length(img)
        i
        imageinfo(i,1)
        slice(imageinfo(i,1),imageinfo(i,2),:,:) = img{i};

    end
    
    Ch1 = double(squeeze(slice(:,WorkChannel,:,:)));%-10020;%-double(min(min(min(slice(:,1,:,:)))));
else
    img = data{1,1}(:,1);
    %zinfo = 1:length(img);
        for i = 1:length(img)
        i
        %imageinfo(i,1)
        slice(i,:,:) = img{i};

    end
    Ch1 = double(squeeze(slice(:,:,:)));
end

img = Ch1;
meta.scale = scale;