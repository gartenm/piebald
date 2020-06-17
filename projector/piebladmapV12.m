function out = piebladmapV12(fn,WorkChannel,invert,P3d,centerX,centerY,centerZ)
% v12 implement index mismatch correction

indcor = 1.383/1.518; %n_cell / n_immersion see Bezlyepkina et al for formula and Park et al for malaria values
%indcor = 1;
anglebin3d = 1;
anglebin = 10;
CylinderP = 1;

%load image
if ~exist('fn','var')
    fn = 'cell1_Airyscan Processing.czi';
end

if ~exist('WorkChannel','var')
    WorkChannel=1;
end
WorkChannelString = strcat('_',num2str(WorkChannel));

if ~exist('invert','var')
    invert=0;
end
if ~exist('P3d','var')
    P3d=1;
end

filetype = fn(end-2:end);

if strcmp(filetype,'tif')
    [Ch1,meta] = OpenTifFile(fn,WorkChannel);
    scale = meta.scale;
elseif strcmp(filetype,'czi')
    [Ch1,meta] = OpenCziZStack(fn,WorkChannel);
    scale = meta.scale;
    scale = scale.*1E6;
end

scale(3) = scale(3)*indcor;

if invert == 1
    Ch1 = invertChannel(Ch1);
end


%background is minimum projection of the stack
Ch1Pz_min = squeeze(min(Ch1,[],1));
for i = 1:size(Ch1,1)
    Ch1bk(i,:,:) = Ch1Pz_min;
    %Ch1bk(i,:,:) = imgaussfilt(squeeze(Ch1(25,:,:)),200);
end


Ch1bkCorr = Ch1-Ch1bk;




% figure(19)
% imagesc(squeeze(Ch1bk(25,:,:)))
%props = regionprops(true(size(A)), A, 'WeightedCentroid')
% figure(20)
% imagesc(squeeze(Ch1(25,:,:)-Ch1bk(25,:,:)))

% tot_mass = sum(sum(sum(Ch1)))
% [zz,ii,jj] = ndgrid(1:size(Ch1,1),1:size(Ch1,2),1:size(Ch1,3));
% Cx = sum(ii(:).*Ch1(:))/tot_mass
% Cy = sum(jj(:).*Ch1(:))/tot_mass
% Cz = sum(zz(:).*Ch1(:))/tot_mass



% figure(11)
% imagesc(Px);
% figure(12)
% imagesc(Py);
% figure(13)
% imagesc(Pz);

% figure(9)
% imagesc(squeeze(Ch1bkCorr(25,:,:)))
%scale = 10000;

if ~exist('centerX','var')
    [centerX, centerY, centerZ] = findCenter(Ch1);
end

% if ~exist('centerX','var')
%     I = int16(Ch1bkCorr);
%     lvl = graythresh(I);
%     BW = imbinarize(I,lvl);
%     maskCh1 = Ch1bkCorr;
%     maskCh1(BW==0)=0;
%     Px = squeeze(sum(maskCh1,2));
%     Py = sum(maskCh1,3);
%     Pz = squeeze(sum(maskCh1,1));
%     A = [1:size(Pz,1)];
%     mass = sum(sum(Pz));
%     centerX = sum(sum(Pz,2).*A')/(mass);
%     centerY = sum(sum(Pz,1).*A)/(mass);
%     A = [1:size(Py,1)];
%     centerZ = sum(sum(Py,2).*A')/(mass);
% end
figure(16)
imagesc(squeeze(Ch1bkCorr(round(centerZ),:,:)))
hold on
plot(centerY,centerX,'ro');
hold off


% figure(14)
% clf;
% hold on
% imagesc(squeeze(Ch1bkCorr(:,round(centerX),:)))
% plot(centerY,centerZ,'ro');
% hold off

% figure(15)
% imagesc(sum(Pz,1));

% VgridX = (([ 1:(size(Ch1bkCorr,2)+1) ] -0.5-(centerX))*scale(1));
% VgridY = (([ 1:(size(Ch1bkCorr,3)+1) ] -0.5-(centerY))*scale(2));
phiSaveName = strcat(fn,'_phi.mat');
if exist(phiSaveName,'file')
    load(phiSaveName)
else
    %tic
    [phi] = makePhi2(Ch1bkCorr,centerX,centerY,scale);
    save(phiSaveName,'phi')
    %timephi = toc
end

if CylinderP==1
map = cylinderPz2(Ch1bkCorr,phi,anglebin);
%map = cylinderPz(Ch1bkCorr,phi,edge,centerX,centerY,1);
figure(99)
imagesc(map)
end



if P3d == 1
    %spherical projection
    angleSaveName = strcat(fn,'_theta.mat');
    newtheta = 1;
    
    if exist(angleSaveName,'file')
        load(angleSaveName)
        
        if sum(savescale == scale) == length(scale)
            newtheta = 0;
        end
    end
    
    if newtheta == 1 %calculate new theta if the scale changed
        disp('new theta is calculated')
        %tic
        [r,theta] = sphericalA2 (Ch1bkCorr,centerX,centerY,centerZ,scale);
        %timetheta = toc
        savescale = scale;
        save(angleSaveName,'r','theta','savescale');
        disp('done')
    end
    
    
    %tic
    [mapS,anglesphi,anglesth] = sphericalP(Ch1bkCorr,r,phi,theta,anglebin3d);
    %time3d = toc
    angles = [(-90+anglebin3d/2):anglebin3d:(90-anglebin3d/2)]*pi/180;
    correction = areaInPfs(scale(3), scale(1), angles);
    %mapSc = mapS'*correction';
    %mapSc = bsxfun(@times, mapS, correction');

%     figure(121);
%     imagesc(mapS);
    
    savename3d = strcat(fn(1:length(fn)-4),WorkChannelString,'_sphere.mat');

    save(savename3d,'anglebin3d','mapS','anglesphi','anglesth');
end


%projection of top and bottom slices
polsize = 1.5; %radius of pol map
[mapT, mapB] = topBottomP(Ch1bkCorr,centerX,centerY,centerZ,polsize,scale);

% figure(122);
% imagesc(mapB);
% figure(123);
% imagesc(mapT);
savename2d = strcat(fn(1:length(fn)-4),WorkChannelString,'_TBmap.mat');
save(savename2d,'mapB','mapT','scale')
savename2d = strcat(fn(1:length(fn)-4),WorkChannelString,'_Tmap.tif');
imwrite(uint16(mapT),savename2d);
savename2d = strcat(fn(1:length(fn)-4),WorkChannelString,'_Bmap.tif');
if ~isempty(mapB)
    imwrite(uint16(mapB),savename2d);
end
%save data

WorkChannelString = strcat('_',num2str(WorkChannel));
if CylinderP==1
savename = strcat(fn(1:length(fn)-4),WorkChannelString,'_CylP.tif');
imwrite(uint16(map),savename);
savename2 = strcat(fn(1:length(fn)-4),WorkChannelString,'_meta.txt');
fileID = fopen(savename2,'w');
%savestring = strcat('zscale (um/px):',{' '}, num2str(scale(3)),'\r\n')
fprintf(fileID, strcat('zscale (um/px):_', num2str(scale(3)),'\r\n'));
fprintf(fileID, strcat('deg/pixel:_', num2str(360/size(map,2))));
fclose(fileID);

%stretch dimension for presentation
multipl = 7;
map2 = NaN(size(map,1)*multipl,size(map,2));
for i = 1:size(map,1)
   for k = 1:multipl
        map2((i-1)*multipl+k,:) = map(i,:);
   end
end
savename = strcat(fn(1:length(fn)-4),WorkChannelString,'_CylP_stretch.tif');
imwrite(uint16(map2),savename);
out.map2 = map2;


savename2 = strcat(fn(1:length(fn)-4),WorkChannelString,'_meta_stretch.txt');
fileID = fopen(savename2,'w');
%savestring = strcat('zscale (�m/px):',{' '}, num2str(scale(3)),'\r\n')
fprintf(fileID, strcat('zscale (um/px):_', num2str(scale(3)/multipl),'\r\n'));
fprintf(fileID, strcat('deg/pixel:_', num2str(360/size(map,2))));
fclose(fileID);
%1um scale bar
scale_width = 10;
scaleStretch = 1/(scale(3)/multipl);
scaleStretch = ones(round(scaleStretch),scale_width)*max(max(map));
savename = strcat(fn(1:length(fn)-4),WorkChannelString,'_CylP_stretch_1umbar.tif');
imwrite(uint16(scaleStretch),savename);
%10 deg scale bar
scaleDeg = 45/(360/size(map,2));
scaleDeg = ones(scale_width, round(scaleDeg))*max(max(map));
savename = strcat(fn(1:length(fn)-4),WorkChannelString,'_CylP_stretch_45degbar.tif');
imwrite(uint16(scaleDeg),savename);
%A = map2(end-size(scaleStretch,1)+1:end,1:scale_width);
map2(end-size(scaleStretch,1)+1:end,1:scale_width) = scaleStretch;
%B = map2(end-scale_width+1:end,1:size(scaleDeg,2))
map2(end-scale_width+1:end,1:size(scaleDeg,2)) = scaleDeg;


savename = strcat(fn(1:length(fn)-4),WorkChannelString,'_CylP_stretch_Wscale.tif');
imwrite(uint16(map2),savename);

out.centerX = centerX;
out.centerY = centerY;
out.centerZ = centerZ;
out.map = map;
out.map3 = map2;
savename = strcat(fn(1:length(fn)-4),WorkChannelString,'.mat');
save(savename,'out');

end

end

function out = cylinderPz(data,cx,cy,wedge)
%3d data as z,x,y
%calculate for each pixle the angle with respect to the center
gridX = [ 1:size(data,2) ] -0.5-(cx);
gridY = [ 1:size(data,3) ] -0.5-(cy);
angleGrid = NaN(length(gridX),length(gridY));
for i = 1: length(gridY)
    for k = 1: length(gridX)
%         angleGrid(i,k) = atan(gridY(i)/gridX(k))*180/pi;
        if gridY(i)> 0
                angleGrid(i,k) = atan(gridX(k)/gridY(i))*180/pi;
        elseif gridY(i) < 0
                if gridX(k) >= 0
                    angleGrid(i,k) = (atan(gridX(k)/gridY(i))+pi)*180/pi;
                else
                    angleGrid(i,k) = (atan(gridX(k)/gridY(i))-pi)*180/pi;
                end
        elseif gridY(i) == 0
                if gridX(k) >= 0
                    angleGrid(i,k) = (pi/2)*180/pi;
                else
                    angleGrid(i,k) = (-pi/2)*180/pi;
                end
        end
            
    end
end
figure(98)
imagesc(angleGrid)
angles = [-180:wedge:180];
map = NaN(size(data,1),length(angles)-1);
for i = 1: size(data,1)
    slice = squeeze(data(i,:,:));
%     figure(101);
%     imagesc(slice)
    for k = 1:length(angles)-1
        

        mask = (angles(k)<=angleGrid) & (angleGrid<angles(k+1));
%         figure(100)
%         imagesc(mask)
        
        map(i,k)=max(max(slice.*mask));
%         figure(99)
%         imagesc(map)
    end
end



out = map;


end

function map = cylinderPz2(data,phi,anglebin)
    edge = phi.edge;
    phiEdge = edge.phi;
    phi = phi.phi;
    anglesphi = gpuArray(([-180:anglebin:180]));
%     phiEdge = phi;
%     h1=edge.data;
%     h2 = phi(1:round(edge.cx),round(edge.cy),2)+360;
%     phiEdge(1:round(edge.cx),round(edge.cy),1) = h1;
%     phiEdge(1:round(edge.cx),round(edge.cy),2) = h2;
%     figure(2);
%     imagesc(squeeze(phiEdge(:,:,1)));
    
    phi = gpuArray((phi));
    phiEdge = gpuArray((phiEdge));
    data = gpuArray((data));
    map = ((zeros(size(data,1),length(anglesphi)-1)));
    %tic
    for k = 1:length(anglesphi)-1
            if anglesphi(k+1)>min(edge.data)
                maskphi = (anglesphi(k)<=phiEdge(:,:,2)) & (phiEdge(:,:,1)<anglesphi(k+1));
            else
                maskphi = (anglesphi(k)<=phi(:,:,2)) & (phi(:,:,1)<anglesphi(k+1));
            end
        for i = 1: size(data,1)
            slice = squeeze(data(i,:,:));
%            mask = maskphi.*maskth;
            maskedData = gather(slice(maskphi==true));
            if ~isempty(maskedData)
                map(i,k)=max(maskedData);
            else
                map(i,k)=0;
            end
        end
    end
    %timepz = toc
    %map = gather(map);
end

function [r,theta] = sphericalA2 (data,cx,cy,cz,scale)
%     cx = round(cx)+0.5;
%     cy = round(cy)+0.5;
%     gridZ = gpuArray(([ 1:(size(data,1)) ]-(cz))*scale(3));
%     gridX = gpuArray(([ 1:(size(data,2)) ]-(cx))*scale(1));
%     gridY = gpuArray(([ 1:(size(data,3)) ]-(cy))*scale(2));
%     
%     
%     VgridZ = gpuArray(([ 1:(size(data,1)+1) ] -0.5-(cz))*scale(3));
%     VgridX = gpuArray(([ 1:(size(data,2)+1) ] -0.5-(cx))*scale(1));
%     VgridY = gpuArray(([ 1:(size(data,3)+1) ] -0.5-(cy))*scale(2));
%     
%     Vphi = gpuArray(NaN(length(VgridX),length(VgridY)));
%     phi = gpuArray(NaN(length(gridX),length(gridY),2));
%     Vr = gpuArray(NaN(length(VgridZ),length(VgridX),length(VgridY)));
%     r = gpuArray(NaN(length(gridZ),length(gridX),length(gridY),2));
%     Vtheta = gpuArray(NaN(length(VgridZ),length(VgridX),length(VgridY)));
%     theta = gpuArray(NaN(length(gridZ),length(gridX),length(gridY),2));
    
    gridZ = (([ 1:(size(data,1)) ]-(cz))*scale(3));
    gridX = (([ 1:(size(data,2)) ]-(cx))*scale(1));
    gridY = (([ 1:(size(data,3)) ]-(cy))*scale(2));
    
    
    VgridZ = (([ 1:(size(data,1)+1) ] -0.5-(cz))*scale(3));
    VgridX = (([ 1:(size(data,2)+1) ] -0.5-(cx))*scale(1));
    VgridY = (([ 1:(size(data,3)+1) ] -0.5-(cy))*scale(2));
    

    Vr = (NaN(length(VgridZ),length(VgridX),length(VgridY)));
    r =(NaN(length(gridZ),length(gridX),length(gridY),2));
    Vtheta = (NaN(length(VgridZ),length(VgridX),length(VgridY)));
    theta = (NaN(length(gridZ),length(gridX),length(gridY),2));
    
    %[phi,edge] = makePhi(VgridX,VgridY,cx,cy);

    
    for xi = 1:length(VgridX)
        for yi = 1:length(VgridY)
            for zi = 1:length(VgridZ)
                Vr(zi,xi,yi) = sqrt(VgridX(xi)^2+VgridY(yi)^2+VgridZ(zi)^2);
            end
        end
    end
    for xi = 1:(length(VgridX)-1)
        for yi = 1:(length(VgridY)-1)
            for zi = 1:(length(VgridZ)-1)
                r(zi,xi,yi,:) = [min(min(min(Vr(zi:(zi+1),xi:(xi+1),yi:(yi+1))))),max(max(max(Vr(zi:(zi+1),xi:(xi+1),yi:(yi+1)))))];
                %r(zi,xi,yi,2) = max(max(max(Vr(zi:(zi+1),xi:(xi+1),yi:(yi+1)))));
            end
        end
    end
    
%     figure(199)
%     imagesc(squeeze(r(round(cz+0),:,:)));
    for xi = 1:length(VgridX)
        for yi = 1:length(VgridY)
            for zi = 1:length(VgridZ)
                Vtheta(zi,xi,yi) = acos(VgridZ(zi)/Vr(zi,xi,yi))*180/pi;
            end
        end
    end
    for xi = 1:(length(VgridX)-1)
        for yi = 1:(length(VgridY)-1)
            for zi = 1:(length(VgridZ)-1)
                theta(zi,xi,yi,:) = [min(min(min(Vtheta(zi:(zi+1),xi:(xi+1),yi:(yi+1))))),max(max(max(Vtheta(zi:(zi+1),xi:(xi+1),yi:(yi+1)))))];
                %theta(zi,xi,yi,2) = max(max(max(Vtheta(zi:(zi+1),xi:(xi+1),yi:(yi+1)))));
            end
        end
    end
    
%     figure(200)
%     imagesc(squeeze(theta(:,round(cx),:,1)));
end
function [map,anglesphi,anglesth] = sphericalP(data,r,phi,theta,anglebin)
   % tic
    edge = phi.edge;
    phiEdge = edge.phi;
    phi = phi.phi;
    
    N = size(data,1)*size(data,2)*size(data,3);
    %anglesphi = gpuArray([-180:anglebin:180]);
    anglesphi = gpuArray([-180:anglebin:180]);
    anglesth = gpuArray([0:anglebin:180]);
    map = (NaN(length(anglesth)-1,length(anglesphi)-1));
    mapk = (NaN(1,length(anglesphi)-1));
    %mask = NaN(size(data));
    %phi = gpuArray(shape(phi,[N,2]));

    theta = gpuArray(reshape(theta,[N,2]));
    
    
%     anglesphi = ([-180:anglebin:180]);
%     anglesth = ([0:anglebin:180]);
%     map = (NaN(length(anglesth)-1,length(anglesphi)-1));
%     %mask = NaN(size(data));
%     phi3d = (NaN(size(data,1),size(data,2),size(data,3),2));
%     phi3dEdge = phi3d;
%     theta = (theta);
%     phi = (phi);
%     data = (data);
    
%     phiEdge = phi;
%     
%     h1=edge.data;
%     h2 = phi(1:round(edge.cx),round(edge.cy),2)+360;
%     phiEdge(1:round(edge.cx),round(edge.cy),1) = h1;
%     phiEdge(1:round(edge.cx),round(edge.cy),2) = h2;
    
    
    
    phi3dH = (NaN(size(data,1),size(data,2),size(data,3)));
    phi3dL = phi3dH;
    phi3dEdgeH = phi3dH;
    phi3dEdgeL = phi3dH;
    for j = 1:size(data,1)
        phi3dL(j,:,:) = phi(:,:,1);
        phi3dH(j,:,:) = phi(:,:,2);
        phi3dEdgeL(j,:,:) = phiEdge(:,:,1);
        phi3dEdgeH(j,:,:) = phiEdge(:,:,2);
    end
%     figure(1002);
%     imagesc(squeeze(phiEdge(:,:,1)-phiEdge(:,:,2)))
    
    
    phi3dL = gpuArray(reshape(phi3dL,[N,1]));
    phi3dH = gpuArray(reshape(phi3dH,[N,1]));
    phi3dEdgeL = gpuArray(reshape(phi3dEdgeL,[N,1]));
    phi3dEdgeH = gpuArray(reshape(phi3dEdgeH,[N,1]));
    data = reshape(data,[N,1]);
%     maskphi4d = false(length(anglesphi)-1,size(data,1),size(data,2),size(data,3));
%     
%     for k = 1:length(anglesphi)-1
%         maskphi4d(k,:,:,:) = (anglesphi(k)<=phi3d) & (phi3d<anglesphi(k+1));
%     end
    maskphi = gpuArray((false(N,1)));
    edgedata = edge.data;

    %maskphi = sparse(maskphi);
    for i = 1: length(anglesth)-1
        %slice = squeeze(data(i,:,:));
        %tic
        maskth = (anglesth(i)<=theta(:,2)) & (theta(:,1)<anglesth(i+1));
        %t1 = toc
%         figure(101);
%         imagesc(squeeze(max(maskth,[],1)))
%         figure(102);
%         imagesc(squeeze(maskth(:,145,:)))
        %tic
        
        for k = 1:length(anglesphi)-1

        if anglesphi(k+1)>min(edgedata)
            maskphi(:) = (anglesphi(k)<=phi3dEdgeH) & (phi3dEdgeL<anglesphi(k+1));
        else
            maskphi(:) = (anglesphi(k)<=phi3dH) & (phi3dL<anglesphi(k+1));
        end
            %tic
%             if anglesphi(k+1)>min(edge.data)
%                 maskphi = (anglesphi(k)<=phi3dEdge(:,:,:,2)) & (phi3dEdge(:,:,:,1)<anglesphi(k+1));
%             else
%                 maskphi = (anglesphi(k)<=phi3d(:,:,:,2)) & (phi3d(:,:,:,1)<anglesphi(k+1));
%             end
            %t2 = toc
%             maskphi2d = (anglesphi(k)<=phi(:,:,2)) & (phi(:,:,1)<anglesphi(k+1));
%             figure(100)
%             imagesc(maskphi2d)
            
%             maskphi = (anglesphi(k)<=phi) & (phi<anglesphi(k+1));
%             for j = 1:size(data,1)
%                 mask(j,:,:) = squeeze(maskth(j,:,:)).*maskphi;
%             end
            %tic
            mask = (maskphi(:) & maskth);
            %t3 = toc
%             mask = squeeze(maskphi4d(k,:,:,:)).*maskth;
%             figure(100)
%             imagesc(maskphi)
            %mask = gather(mask);
          %tic
            %idx = (mask==true);
            idx = find(mask);
          % t4 = toc
            maskedData = (data(idx));
            %maskedData = nonzeros(data.*mask);
           
           %maskedData = gather(maskedData);
           
%             tic
%             maskedData2 = maskedData(maskedData>0);
%                         toc
            %MaxData = max(maskedData);
           % tic
            if ~isempty(maskedData)
                mapk(k)=(max(maskedData));
            else
                mapk(k)=0;
            end
            %t5 = toc

    %         figure(99)
    %         imagesc(map)
        end
        map(i,:) = mapk;
        %tparfor = toc
    end
    
    figure(104)
    imagesc(map)
    %map = gather(map);
   % toc
end

function [phis] = makePhi2(data,cx,cy,scale)
    gridX = (([ 1:(size(data,2)) ]-(cx))*scale(1));
    gridY = (([ 1:(size(data,3)) ]-(cy))*scale(2));
    VgridX = (([ 1:(size(data,2)+1) ] -0.5-(cx))*scale(1));
    VgridY = (([ 1:(size(data,3)+1) ] -0.5-(cy))*scale(2));
    
    Vphi = (NaN(length(VgridX),length(VgridY)));
    phi = (NaN(length(gridX),length(gridY),2));
    for i = 1:length(VgridX)
        Vxi = VgridX(i);
        for j = 1:length(VgridY)
            Vphi(i,j) = atan2(VgridY(j),Vxi)*180/pi;
        end
    end
    
    for i = 1:(length(VgridX)-1)
        for j = 1:(length(VgridY)-1)
            phi(i,j,:) = [min(min(Vphi(i:(i+1),j:(j+1)))),max(max(Vphi(i:(i+1),j:(j+1))))];
            %phi(i,j,2) = max(max(Vphi(i:(i+1),j:(j+1))));
        end
    end
    
    phiEdge = phi;
    h=phi(1:round(cx),round(cy),2)-360;
    edge.data = phi(1:round(cx),round(cy),2);
    edge.edgeHigh = phi(1:round(cx),round(cy),1)+360;
    phiEdge(1:round(cx),round(cy),2) = edge.edgeHigh;
    
    edge.edgeLow= phi(1:round(cx),round(cy),2);
    phiEdge(1:round(cx),round(cy),1) = edge.edgeLow;


    edge.data1New = h;
    edge.data2New = phi(1:round(cx),round(cy),1);
    
    
    edge.cx = cx;
    edge.cy = cy;
    phis.edge = edge;
    
    phi(1:round(cx),round(cy),2)=phi(1:round(cx),round(cy),1);
    phi(1:round(cx),round(cy),1) = h;
    figure(198)
    imagesc(squeeze(phi(:,:,2)-phi(:,:,1)))
    
    phis.phi = phi;
    phis.edge.phi = phiEdge;
    
end    
    
function correction = areaInPfs(zres, xyres, angles)
    r = (zres*xyres)./sqrt(((xyres.*cos(angles)).^2)+((zres.*sin(angles)).^2));
    A = xyres.*r.*pi;
    correction = max(A)./A;

end

function [mapT, mapB] = topBottomP(data,cx,cy,cz,polsize,scale)

    Zoffset = 10;
    sx = (polsize/scale(1));
    %sxm = sx/2;
    sy = (polsize/scale(2));
    %sym = sy/2;
%     x1 = round(cx-sx);
%     x2 = round(cx+sx);
%     y1 = round(cy-sy);
%     y2 = round(cy+sy);
    x1 = 1;
    x2 = size(data,2);
    y1 = 1;
    y2 = size(data,3);

%     mapB = squeeze(max(data(1:(cz-1),:,:),[],1));
%     mapT = squeeze(max(data((cz+1):end,:,:),[],1));
    
    mapB = squeeze(max(data(1:round(cz-Zoffset),x1:x2,y1:y2),[],1));
    mapT = squeeze(max(data(round(cz+Zoffset):end,x1:x2,y1:y2),[],1));
%     rmax = sqrt((sx)^2+(sy)^2);
%     %polmask = NaN(sx,sy);
%     for xi=1:size(mapB,1)
%         for yi=1:size(mapB,2)
%             r = sqrt((xi-sx)^2+(yi-sy)^2);
%             if r>rmax
%                 mapB(xi,yi)=NaN;
%                 mapT(xi,yi)=NaN;
%             end
%         end
%     end
end

function out = invertChannel(Ch)
    out = -Ch+6553;
end
