function simplyFitHistogram

fnlist = dir('*PVMPPMdist.mat');

binwidth = 0.5;
binedges = 0:binwidth:60;
centers = binedges(1:(end-1))+binwidth/2;
N = 1;
d1 = zeros(N,1);
d2 = d1;
s1 = d1;
s2 = d1;
a = d1;
medDist = d1;
sseres = d1;
datastr = {};
for i = 1:length(fnlist)
    datastr{i} = load(fnlist(i).name);
end

for Nstrp=1:N
    A = centers*0;
    B = [];
    Nstrp
    idx = [1:length(fnlist)];
    %randi([1 length(fnlist)],[1 length(fnlist)]);

    for i = idx
        %load(fnlist(i).name);
        %fnlist(i).name;
        %max(Ch1Ch2Dist(:,6));
        dataset = datastr{1, i}.Ch1Ch2Dist(:,6);
        scale = datastr{1, i}.scale;
        [N1,edges] = histcounts(dataset,binedges);
        %[N1,edges] = histcounts(Ch1Ch2Dist(:,6),binedges);
        N = N1*scale;
        A = A+N;
        %B = [B; Ch1Ch2Dist(:,6)];
        B = [B; dataset];
    end

    Avector = [];

    for i=1:length(A)
        segmentL = A(i); %hight of the scaled bar in nm
        dist = centers(i); % value of the membrane distance
        for k=1:round(segmentL)
            Avector = [Avector dist];
        end
    end
%     logcenters = log(centers);
%     PPMPVMnormcdf = cumsum(A);
%     PPMPVMnormcdf = PPMPVMnormcdf/max(PPMPVMnormcdf);

    [fitresult, gof] = AnalyzeDistro(Avector);
    
    d1(Nstrp) = exp(fitresult.mu1);
    s1(Nstrp) = exp(fitresult.sigma1);
    d2(Nstrp) = exp(fitresult.mu2);
    s2(Nstrp) = exp(fitresult.sigma2);
    a(Nstrp) = fitresult.a;
    medDist(Nstrp) = median(Avector);
    sseres(Nstrp) = gof.sse;
end

ci = confint(fitresult)
val = coeffvalues(fitresult)
expci = exp(ci)
expval = exp(val)

sumA = sum(A)/1000

% pdd1 = fitdist(d1,'Normal')
% pdd2 = fitdist(d2,'Normal')
% pda = fitdist(a,'Normal')
% pdmedDist = fitdist(medDist,'Normal')
% d1ci = paramci(fitdist(d1,'Normal'))
% d2ci = paramci(fitdist(d2,'Normal'))
% aci = paramci(fitdist(a,'Normal'))


% d1result = CI95(d1)
% d2result = CI95(d2)
% s1result = CI95(s1)
% s2result = CI95(s2)
% aresult = CI95(a)
% medresult = CI95(medDist)


save('simplefit.mat','fitresult','gof','ci','expci','val','expval','sumA')

% figure(324)
% clf;
% hold on;
% title('mean 1')
% hist(d1);
% 
% 
% %hist(d2);
% figure(325);
% clf;
% hold on;
% title('mean 2')
% hist(d2);
% 
% figure(326);
% clf;
% hold on;
% title('relative contribution')
% hist(a)
% 
% figure(327);
% clf;
% hold on;
% title('median')
% hist(medDist)