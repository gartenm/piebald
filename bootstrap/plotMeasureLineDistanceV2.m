function plotMeasureLineDistanceV2 %[A,centers,Ng,histdata,B,gold] = 

fnlist = dir('*PVMPPMdist.mat');

binwidth = 2;
binedges = 0:binwidth:70;
centers = binedges(1:(end-1))+binwidth/2;
A = centers*0;
B = [];
for i = 1:length(fnlist)
    load(fnlist(i).name);
    fnlist(i).name
    max(Ch1Ch2Dist(:,6))
    [N1,edges] = histcounts(Ch1Ch2Dist(:,6),binedges);
    N = N1*scale;
    A = A+N;
    B = [B; Ch1Ch2Dist(:,6)];
    
end

%A(A>100) = [];

figure(237);
%[counts,centers] = hist(A,30);
bar(centers,A);
%hist(A,30);
ylabel('length of PV (nm)');
xlabel('PVM - PPM distance (nm)');
areaA = sum(binwidth*A);


save('PVMPPMdist_result.mat','A','centers')

figure(238);
ecdf(B);
ylabel('cummulative counts');
xlabel('PVM - PPM distance (nm)')
%Ax = A/27;


%cut the histogram in frequency data
    Avector = [];
    
    for i=1:length(A)
        segmentL = A(i); %hight of the scaled bar in nm
        dist = centers(i); % value of the membrane distance
        for k=1:round(segmentL)
            Avector = [Avector dist];
        end
    end
    logcenters = log(centers);
    PPMPVMnormcdf = cumsum(A);
    PPMPVMnormcdf = PPMPVMnormcdf/max(PPMPVMnormcdf);
figure(260);
clf;
hold on;
h1 = histogram(Avector)
h1.BinWidth = binwidth;
ylabel('length of PV (nm)');
xlabel('PVM - PPM distance (nm)');
%golddistFN = 'goldDist.mat';
%golddistFN = 'gold2.mat';
golddistFN = 'fakegold.mat';
if isfile(golddistFN)

    load(golddistFN);
    figure(239);
    [Ng,edges] = histcounts(gold,binedges);
    bar(centers,Ng);
    areagold = sum(binwidth*Ng);
    %hist(A,30);
    ylabel('number gold');
    xlabel('PVM - PPM distance (nm)')
    
    figure(240);
    clf;
    hold on;
    
    histdata =  [(A/max(A)); (Ng/max(Ng))]';
    %histdata =  [(A/areaA); (Ng/areagold)]'
    %X = [centers];
    bar(centers,histdata);
    %bar(centers,Ng/areagold,'o');
    ylabel('probability density');
    xlabel('PVM - PPM distance (nm)');
    

    
    [goldnormcdf, goldx] = ecdf(gold);
    
    figure(8);
    clf;
    hold on;
    plot(centers,PPMPVMnormcdf);
    plot(goldx,goldnormcdf);
    

    
    ecdf(Avector);
    ecdf(gold);
    
    
    logcenters = log(centers);
    loggoldx = log(goldx);
    figure(9);
    clf;
    hold on;
    plot(logcenters,PPMPVMnormcdf);
    plot(loggoldx,goldnormcdf);
    
    save('distroVars.mat','Ng','gold','A','B','Avector','centers','logcenters','PPMPVMnormcdf','goldx','loggoldx','goldnormcdf')
    
    %cCount = centers(centers<20);
    TotLclose = sum(A(centers<=20))/1000;
    TotLfar = sum(A(centers>20))/1000;
    NgoldClose = sum(gold<=20);
    NgoldFar = sum(gold>20);
    
    ratioFar = NgoldFar/TotLfar
    ratioClose = NgoldClose/TotLclose
    
    median(gold)
    %for i=1:length(cCount)
        
    
%    normcdf
    
    kstest2(Avector,gold)
else 
    
     save('distroVars.mat','A','B','Avector','centers','logcenters','PPMPVMnormcdf')
    
end

