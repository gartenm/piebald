function bootstrapallPVMPPMdist

%golddistFN = 'goldDist.mat';
golddistFN = 'alldara01resample.mat';
%golddistFN = 'fakegold.mat';
load(golddistFN);
gold = Avector;
N = 10000;
d1 = zeros(N,1);
d2 = d1;
a = d1;
medDist = d1;
for Nstrp=1:N
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
aresult = CI95(a)
medresult = CI95(medDist)


save('bootstrapPVMPPMpointbypoint10k.mat','d1result','d2result','aresult','medresult')

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