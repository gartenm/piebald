function bootstrapallPVMPPMdistFlexiblePVlength

%golddistFN = 'goldDist.mat';
golddistFN = 'alldara01resample.mat';
%golddistFN = 'fakegold.mat';
load(golddistFN);

load('distroVars.mat', 'A')
load('distroVars.mat', 'centers')

Avector = [];
Arc = (18+9)*2; %length of membrane segments corresponding to a possible binding site of a gold particle

for i=1:length(A)
    segmentL = A(i)/Arc; %hight of the scaled bar in nm
    dist = centers(i); % value of the membrane distance
    for k=1:round(segmentL)
        Avector = [Avector dist];
    end
end
%     logcenters = log(centers);
%     PPMPVMnormcdf = cumsum(A);
%     PPMPVMnormcdf = PPMPVMnormcdf/max(PPMPVMnormcdf);

    length(Avector)
    
gold = Avector;
N = 10000;
d1 = zeros(N,1);
d2 = d1;
s1 = d1;
s2 = d1;

a = d1;
medDist = d1;
parfor Nstrp=1:N
    Nstrp
    idx = randi([1 length(gold)],[1 length(gold)]);

    %for i = idx
        
        bgold = gold(idx);
%         figure(239);
%         hist(bgold);
%         figure(240);
%         ecdf(bgold);
%         
%         [Ng,edges] = histcounts(gold,binedges);
%         bar(centers,Ng);
        %areagold = sum(binwidth*Ng);
        %hist(A,30);
%         ylabel('number gold');
%         xlabel('PVM - PPM distance (nm)')
% 
%         figure(240);
%         clf;
%         hold on;
    [fitresult, gof] = AnalyzeDistro(bgold);

    d1(Nstrp) = exp(fitresult.mu1);
    d2(Nstrp) = exp(fitresult.mu2);
    s1(Nstrp) = exp(fitresult.sigma1);
    s2(Nstrp) = exp(fitresult.sigma2);
    a(Nstrp) = fitresult.a;
    medDist(Nstrp) = median(bgold);
end


% pdd1 = fitdist(d1,'Normal')
% pdd2 = fitdist(d2,'Normal')
% pda = fitdist(a,'Normal')
% pdmedDist = fitdist(medDist,'Normal')
% d1ci = paramci(fitdist(d1,'Normal'))
% d2ci = paramci(fitdist(d2,'Normal'))
% aci = paramci(fitdist(a,'Normal'))

d1result = CI95(d1)
d2result = CI95(d2)
s1result = CI95(s1)
s2result = CI95(s2)
aresult = CI95(a)
medresult = CI95(medDist)


save('bootstrapPVMPPMArc10kNew.mat','d1result','d2result','s1result','s2result','aresult','medresult')

figure(324)
clf;
hold on;
title('mean 1')
hist(d1);


%hist(d2);
figure(325);
clf;
hold on;
title('mean 2')
hist(d2);

figure(326);
clf;
hold on;
title('relative contribution')
hist(a)

figure(327);
clf;
hold on;
title('median')
hist(medDist)