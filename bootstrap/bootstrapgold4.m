function bootstrapgold4
% only a is varied for the fit
%golddistFN = 'goldDist.mat';
golddistFN = 'gold2.mat';
%golddistFN = 'fakegold.mat';
load(golddistFN);

N = 10000;
d1 = zeros(N,1);
d2 = d1;
a = d1;
medDist = d1;
sse = d1;
aic = d1;
bic = d1;
np = 3;

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
    [fitresult, gof] = AnalyzeDistroGold3(bgold)

%    d1(Nstrp) = exp(fitresult.mu1);
    %d2(Nstrp) = exp(fitresult.mu2);
   
    a(Nstrp) = fitresult.a;
    medDist(Nstrp) = median(bgold);
    sse(Nstrp)=gof.sse;
    aic(Nstrp)=2*np-2*log(gof.sse);
    bic(Nstrp)=np*log(length(gold))-2*log(gof.sse);
end

Rsse.sse = sse;
Rsse.aic = aic;
Rsse.bic = bic;
%pdd1 = fitdist(d1,'Normal')
%pdd2 = fitdist(d2,'Normal')
pda = fitdist(a,'Normal')
pdmedDist = fitdist(medDist,'Normal')

%d1result = CI95(d1)
%d2result = CI95(d2)
aresult = CI95(a)
medresult = CI95(medDist)


save('goldbootstrap410ksseNew.mat','aresult','medresult','Rsse')

% d1ci = paramci(fitdist(d1,'Normal'))
% d2ci = paramci(fitdist(d2,'Normal'))
% aci = paramci(fitdist(a,'Normal'))


%save('goldbootstrap.mat','pdd1','pda','pdmedDist')

figure(324)
clf;
hold on;
title('mean 1')
%hist(d1);


%hist(d2);
figure(325);
clf;
hold on;
title('mean 2')
%hist(d2);

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

figure(328);
clf;
hold on;
title('sse')
hist(Rsse.sse)

figure(329);
clf;
hold on;
title('aic')
hist(Rsse.aic)

figure(330);
clf;
hold on;
title('bic')
hist(Rsse.bic)