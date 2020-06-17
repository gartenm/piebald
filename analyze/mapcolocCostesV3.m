function cv1v2 = mapcolocCostesV3
%works with spherical mollweide projections
% fn1 = '20 cells_Airyscan Processing-Scene-06_1_sphere.mat_mollweideP.tif';
% fn2 = '20 cells_Airyscan Processing-Scene-06_2_sphere.mat_mollweideP.tif';

%fn1list(1).name = '16 cells_Airyscan Processing-Scene-02_1_sphere.mat_mollweideP.tif';
% fn2list(1).name = '20 cells_Airyscan Processing-Scene-06_2_sphere.mat_mollweideP.tif';

fn1list = dir('*1_sphere.mat_mollweideP.tif'); 
% fn2list = dir('*_2.mat');
figure(274);
clf;
hold on;
xbins = [-1:0.1:1];
cv1v2res = [];
rcostesResult = [];
mandersResult = [];
rcostesResultTotal = [];
mandersResultTotal = [];
for iFile = 1: length(fn1list)
    fn1 = fn1list(iFile).name
    basename = fn1(1:(length(fn1)-27));
    fn2 = strcat(basename,'2_sphere.mat_mollweideP.tif');

    %load(fn1);
    map1o = double(imread(fn1));

    %load(fn2);
    map2o = double(imread(fn2));
    %cz = round(out.centerZ);
    %238:270
    %17:22
    idx0 = (map1o == 0)&(map2o == 0);
    map1 = map1o;
    map2 = map2o;
    map1(idx0)=NaN;
    map2(idx0)=NaN;
%     map1(isnan(map1))=20000;
%     map2(isnan(map2))=20000;
%     SnanMap1 = sum(sum(isnan(map1)));
%     SnanMap2 = sum(sum(isnan(map2)));
                    
    savename = strcat(fn1,'_coloc.mat');
    if exist(savename,'file')
        reload = load(savename);
        mask = reload.out.mask;
    else
        mask = maskDefects(map1,map2);
    end
    map1(mask==1)=NaN;
    map2(mask==1)=NaN;

    figure(232);
    imagesc(map1)
    figure(233);
    imagesc(map2)
    
    

    N = size(map1,1)*size(map1,2);
    C1 = reshape(map1,[N,1]);
    C2 = reshape(map2,[N,1]);
    C1(isnan(C1))=[];
    C2(isnan(C2))=[];
    [rTotal,M12Total] = plotPearsonCorr(C1,C2,0,0);
    [lvl1,lvl2] = costes(C1,C2);
    [rcostes,M12costes] = plotPearsonCorr(C1,C2,lvl1,lvl2);
    
%     f = 1;
%     medC1L = mean(C1(C1<lvl1))*f;
%     medC1H = mean(C1(C1>lvl1))*f;
%     map1sc = ((map1-medC1L)./(medC1H-medC1L).*2)-1;
% 
%     medC2L = mean(C2(C2<lvl2))*f;
%     medC2H = mean(C2(C2>lvl2))*f;
%     map2sc = ((map2-medC2L)./(medC2H-medC2L).*2)-1;
%     
%     map1sc(map1sc > 1) = 1;
%     map1sc(map1sc < -1) = -1;
%     map2sc(map2sc > 1) = 1;
%     map2sc(map2sc < -1) = -1;


    
    map1m = map1 - min(min(map1));
    map2m = map2 - min(min(map2));
    
    
%     map1sc = map1m./max(max(map1m));
%     map2sc = map2m./max(max(map2m));
    
    map1sc = (map1m./max(max(map1m))).*2-1;
    map2sc = (map2m./max(max(map2m))).*2-1;
% %   
%     map1sc(map1 < lvl1) = -1;
%     map2sc(map2 < lvl2) = -1;
       
    
    
%     lvl1 = graythresh(map1sc);
%     lvl2 = graythresh(map2sc);
    
%Yn = size(Y,1)*size(Y,2);



    
    
%     lvl1m = multithresh(map1sc,2)
%     lvl2m = multithresh(map2sc,2)
%     lvl1 = lvl1m(1);
%     lvl2 = lvl2m(1);
    
    mapbin1 = (map1 > lvl1);
    figure(244);
    imagesc(mapbin1)
    mapbin2 = (map2 > lvl2);
    figure(245);
    imagesc(mapbin2)
    
    mapbin = (map1 > lvl1) | (map2 > lvl2);
    figure(243);
    imagesc(mapbin)
    
    mapbinZhist = sum(mapbin,2);
%     istart = find(mapbinZhist>0,1,'first')+4;
%     iend = find(mapbinZhist>0,1,'last')-12;
    cz = size(map1,1)/2;
    istart = cz-6;
    iend = cz+6;
    cv1v2 = [];

    figure(252);
    imagesc(map1sc)
    figure(253);
    imagesc(map2sc)
    
    
    for i = istart:iend%25:size(map1,1)
%         v1 = (map1sc(i,:)-lvl1).*(lvl1m(2)-lvl1m(1));
%         v2 = (map2sc(i,:)-lvl2).*(lvl2m(2)-lvl2m(1));
        v1 = (map1sc(i,:));
        v2 = (map2sc(i,:));
        v1 = v1(~isnan(v1));
        v2 = v2(~isnan(v2));
%         v1 = (map1sc(i,:));%-lvl1;
%         v2 = (map2sc(i,:));%-lvl2;
%         v1 = (map1sc(i,:).*2)-1;%-lvl1;
%         v2 = (map2sc(i,:).*2)-1;%-lvl2;
        v1norm = v1./norm(v1);
        v2norm = v2./norm(v2);
        cv1v2 = [cv1v2 dot(v1norm,v2norm)];
        %cv1v1 = dot(v1norm,v1norm);
        %cv2v2 = dot(v2norm,v2norm);
    end
    figure(234);
    plot(cv1v2);
    figure(274);
    hist(cv1v2,xbins);
    cv1v2res = [cv1v2res, cv1v2];
    rcostesResult = [rcostesResult, rcostes];
    mandersResult = [mandersResult, M12costes];
    rcostesResultTotal = [rcostesResultTotal, rTotal];
    mandersResultTotal = [mandersResultTotal, M12Total];
    %savename = strcat(fn1,'coloc.mat');
    out.mask = mask;
    out.map1 = map1;
    out.map2 = map2;
    out.rcostes = rcostes;
    out.manders = M12costes;
    out.rcostesTotal = rTotal;
    out.mandersTotal = M12costes;
    out.cv1v2 = cv1v2;
    out.costesT1 = lvl1;
    out.costesT2 = lvl2;
    save(savename,'out')
    
    text_str{1} = strcat('R: ', num2str(rTotal));
    text_str{2} = strcat('Manders: ', num2str(M12Total));
    position = [20 20; 20 50];
    tx = size(map1,1)/2;
    ty = size(map1,2)/2
    I = zeros(tx,ty);
    
    im = addtext(I, text_str,position);
    
    %RGB = insertText(I,position,text_str);

    infolayer = zeros(size(map1,1),size(map1,2));
    infolayer((tx+1):end,1:ty)= im;
    figure(2131);
    imshow(infolayer);
    
    savename= strcat(basename,'_mollweide_info.tif');
    imwrite(uint16(map1o), savename, 'Compression','none');
    imwrite(uint16(map2o), savename, 'WriteMode', 'append',  'Compression','none');
    imwrite(uint16(infolayer), savename, 'WriteMode', 'append',  'Compression','none');

    
end
    figure(278);
    [counts,centers] = hist(cv1v2res,xbins);
    hist(cv1v2res,xbins);
    ylabel('lines')
    xlabel('correlation')
    
    z = .5.*(log(1+centers) - log(1-centers));
    
    figure(279);
    plot(z,counts,'sk');
    
    ylabel('lines')
    xlabel('fisher z')
    
    se = 1/sqrt(length(v1norm))
    save('mapcolocResult.mat','cv1v2res','rcostesResult','mandersResult','rcostesResultTotal','mandersResultTotal')
    
    figure(214);
    hist(rcostesResult)
    title('Pearson after Costes');
    figure(215);
    hist(mandersResult)
    title('Manders after Costes');
    figure(216);
    hist(rcostesResultTotal)
    title('Pearson no threshold');
    figure(217);
    hist(mandersResultTotal)
    title('Manders no threshold');
    

    
    
    
function [rcostes,M12] = plotPearsonCorr(C1,C2,T1,T2)
    idx = ((C1>T1)|(C2>T2));


    figure(533);
    clf;
    hold on;
    plot(C1(idx),C2(idx),'.');
    xlabel('Channel 1 intensity');
    ylabel('Channel 2 intensity');


    rcostes = pearsonCorr(C1(idx),C2(idx))

    ColocArray1 = C1(idx);
    ColocArray2 = C2(idx);


    [n,c] = hist3([ColocArray1,ColocArray2],[100,100]);
    max(ColocArray1)
    max(ColocArray2)



    % figure(100);
    % imagesc(maskedIM(:,:,imidx,1));
    % figure(102);
    % imagesc(maskedIM(:,:,imidx,2));

    M12 = manders(ColocArray1,ColocArray2)

    figure(101);
    h = pcolor(c{1,1},c{1,2},n');
    set(h, 'edgecolor','none');
    %ylim([0 1]);
    %xlim([0 1]);
    xlabel('channel 1 intensity');
    ylabel('channel 2 intensity');

function mandersCoeff = manders(R,G)
    mandersCoeff = nansum(R.*G)/sqrt(nansum(R.*R)*nansum(G.*G));
    
function mask = maskDefects(C1,C2)
    RGB = zeros(size(C1,1),size(C1,2),3);
    RGB(:,:,1) = C1/(max(max(C1)));
    RGB(:,:,2) = C2/(max(max(C2)));
    fig = figure(891);
    image(RGB)
    title('chose masking area by clicking corners rois, enter when done')
    mask = zeros(size(C1,1),size(C1,2));
    button = 1;
    xprev = [];
    while ~isempty(button)
        %rect = getrect(fig) 
        [x,y,button] = ginput(1)
        if isempty(xprev)
            xprev = x;
            yprev = y;
        else
            mask(min(yprev,y):max(yprev,y),min(xprev,x):max(xprev,x)) = 1;
            xprev = [];
            RGB(:,:,3) = mask;
            image(RGB);
        end
        
    end
    
function im = addtext(im, txt,pos)
% Read example image 
%load('clown.mat'); 
%im = uint8(255*ind2rgb(X,map));
%%Create the text mask 
% Make an image the same size and put text in it 
hf = figure('color','white','units','normalized','position',[.1 .1 .8 .8]); 
image(ones(size(im))); 
set(gca,'units','pixels','position',[5 5 size(im,2) size(im,1)],'visible','off')
% Text at arbitrary position
for i = 1:length(txt)
text('units','pixels','position',pos(:,i),'fontsize',30,'string',txt{i});
end
% Capture the text image 
% Note that the size will have changed by about 1 pixel 
tim = getframe(gca); 
close(hf) 
% Extract the cdata
tim2 = sum(tim.cdata,3);
% Make a mask with the negative of the text 
tmask = tim2==0; 
% Place white text 
% Replace mask pixels with UINT8 max 
im(tmask) = 255; 
%image(im);
%axis off
    