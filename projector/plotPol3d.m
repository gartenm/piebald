function [mapNorth,mapSouth, mapTif] = plotPol3d(fn)

if ~exist('fn','var')
    fn = 'cell1_1_sphere.mat';
end

load(fn);

projectionsize = 1024;
midpoint = projectionsize/2;

minlatitude = 40;

r = midpoint/cos(minlatitude*pi/180);
phi = zeros(projectionsize,projectionsize);
theta = zeros(projectionsize,projectionsize);
    for i = 1:projectionsize
        for j = 1:projectionsize
            phi(i,j) = atan2(j-midpoint,i-midpoint)*180/pi;
            d = sqrt((j-midpoint)^2+(i-midpoint)^2);
            theta(i,j) = real(acos(d/r)*180/pi);
        end
    end
figure(425);
imagesc(phi)
figure(426);
imagesc(theta)

thetaMap = anglesth-90;
ithetaMapN = find(thetaMap>=minlatitude);
ithetaMapS = find(thetaMap<=-minlatitude);


mapNorth = polprojector(thetaMap,ithetaMapN,theta,anglesphi,phi,mapS,projectionsize);
figure(427);
imagesc(mapNorth); 
mapSouth = polprojector(thetaMap,ithetaMapS,theta,anglesphi,phi,mapS,projectionsize);
figure(428);
imagesc(mapSouth); 

% if exist('MollweideLatLong.mat','file')
%     load('MollweideLatLong.mat');
% else
    [pmapLat,pmapLong] = mollweide(projectionsize);
%     save('MollweideLatLong.mat','pmapLat','pmapLong');
% end

map = MapP(pmapLat,pmapLong,mapS,thetaMap,theta,anglesphi,phi);

figure(429);
imagesc(map); 

mapTif = map;
mapTif (isnan(map)) = 0;

savename = strcat(fn,'_mollweideP.tif');
imwrite(uint16(mapTif),savename);
savename = strcat(fn,'_projections.mat');
save(savename,'mapNorth','mapSouth','map')


% for i = 1:(length(ithetaMapN)-1)
%     thSt = thetaMap(ithetaMapN(i));
%     thEnd = thetaMap(ithetaMapN(i+1));
%     thmask = ((theta>thSt)&(theta<=thEnd));
%     for j=1:(length(anglesphi)-1)
%         phiSt = anglesphi(j);
%         phiEnd = anglesphi(j+1);
%         phimask = ((phi>phiSt)&(phi<=phiEnd));
%         mask = thmask & phimask;
% %     figure(428)
% %     imagesc(mask)
%         mapNorth(mask == 1) = mapS(ithetaMapN(i),j);
%     end
%     figure(427)
%     imagesc(mapNorth) 
% end

function map = polprojector(thetaMap,ith,theta,anglesphi,phi,mapS,projectionsize)

map = NaN(projectionsize,projectionsize);
phi = gpuArray(phi);
theta = gpuArray(theta);
for i = 1:(length(ith)-1)
    thSt = abs(thetaMap(ith(i)));
    thEnd = abs(thetaMap(ith(i+1)));
    thmask = ((theta>min([thSt,thEnd]))&(theta<=max([thSt,thEnd])));
    for j=1:(length(anglesphi)-1)
        phiSt = anglesphi(j);
        phiEnd = anglesphi(j+1);
        phimask = ((phi>phiSt)&(phi<=phiEnd));
        mask = thmask & phimask;
%     figure(428)
%     imagesc(mask)
        idx = find(mask);
        map(idx) = mapS(ith(i),j);
    end

end

function [pmapLat,pmapLong] = mollweide(ps)

R = ps/(2*sqrt(2));

pmapLat = NaN(ps,2*ps);
pmapLong = NaN(ps,2*ps);
figure(12);
clf;
hold on;

for yi = 1:ps
    y = ((yi-(ps/2))/(ps/2))*(R*sqrt(2));
    th = asin(y/(R*sqrt(2)));
    plot(yi,th,'ks')
   % asin((2*th+sin(2*th))/pi)
   
    pmapLat(yi,:) = real(asin((2*th+sin(2*th))/pi)*180/pi);
    for xi = 1:(2*ps)
   %     (pi*xi)/(2*R*sqrt(2)*cos(th))
        x = ((xi-ps)/ps)*(2*R*sqrt(2));
        L = real((pi*x)/(2*R*sqrt(2)*cos(th))*180/pi);
        if (L>-180) && (L<180)
            pmapLong(yi,xi) = real((pi*x)/(2*R*sqrt(2)*cos(th))*180/pi);
        end
    end
end

figure(323)
imagesc(pmapLat);
figure(322)
imagesc(pmapLong);

function map = MapP(pmapLat,pmapLong,mapS,thetaMap,theta,anglesphi,phi)

map = NaN(size(pmapLat));
GpmapLat = gpuArray(pmapLat);
GpmapLong = gpuArray(pmapLong);

for i = 1:(length(thetaMap)-1)
    thSt = (thetaMap(i));
    thEnd = (thetaMap(i+1));
    
    thmask = ((GpmapLat>thSt)&(GpmapLat<=thEnd));
    for j=1:(length(anglesphi)-1)
        phiSt = anglesphi(j);
        phiEnd = anglesphi(j+1);
        phimask = ((GpmapLong>phiSt)&(GpmapLong<=phiEnd));
        mask = thmask & phimask;
%     figure(428)
%     imagesc(mask)
        map(mask == 1) = mapS(i,j);
    end

end

function tn1 = thetaNp1(t,p)

tn1 = t-(2*t+sin(2*t)-pi*sin(p))/(2+2*cos(2*t));


