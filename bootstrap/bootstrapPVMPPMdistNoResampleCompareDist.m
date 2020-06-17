function bootstrapPVMPPMdistNoResampleCompareDist
%randomly picks samples without reusing, makes two distributions per
%bootstap

n1 = 34;
n2 = 46;

fnlist = dir('*PVMPPMdist.mat');

%binwidth = 0.5;
%binedges = 0:binwidth:60;
%centers = binedges(1:(end-1))+binwidth/2;
N = 10000;
d11 = zeros(N,1);
d12 = d11;
a1 = d11;
medDist1 = d11;
d21 = zeros(N,1);
d22 = d11;
a2 = d11;
medDist2 = d11;
datastr = {};
n = length(fnlist);


for i = 1:n
    datastr{i} = load(fnlist(i).name);
end

parfor Nstrp=1:N
    %A = centers*0;
    B1 = [];
    B2 = [];
    %Nstrp
    p = randi([1 n],[1 n]);
    %p = randperm(n);
    idx1 = p(1:n1);
    idx2 = p((n1+1):(n2+n1));
    
    for i = idx1
        %load(fnlist(i).name);
        %fnlist(i).name;
        dataset = datastr{1, i}.Ch1Ch2Dist(:,6);
        %max(Ch1Ch2Dist(:,6));
        %[N1,edges] = histcounts(Ch1Ch2Dist(:,6),binedges);
        %N = N1*scale;
        %A = A+N;
        B1 = [B1; dataset];
    end
    for i = idx2
        %load(fnlist(i).name);
        %fnlist(i).name;
        dataset = datastr{1, i}.Ch1Ch2Dist(:,6);
        %max(Ch1Ch2Dist(:,6));
        %[N1,edges] = histcounts(Ch1Ch2Dist(:,6),binedges);
        %N = N1*scale;
        %A = A+N;
        B2 = [B2; dataset];
    end
    %Avector = [];

%     for i=1:length(A)
%         segmentL = A(i); %hight of the scaled bar in nm
%         dist = centers(i); % value of the membrane distance
%         for k=1:round(segmentL)
%             Avector = [Avector dist];
%         end
%     end
%     logcenters = log(centers);
%     PPMPVMnormcdf = cumsum(A);
%     PPMPVMnormcdf = PPMPVMnormcdf/max(PPMPVMnormcdf);

    [fitresult, gof] = AnalyzeDistro(B1);

    d11(Nstrp) = exp(fitresult.mu1);
    d12(Nstrp) = exp(fitresult.mu2);
    a1(Nstrp) = fitresult.a;
    medDist1(Nstrp) = median(B1);
    
    [fitresult, gof] = AnalyzeDistro(B2);

    d21(Nstrp) = exp(fitresult.mu1);
    d22(Nstrp) = exp(fitresult.mu2);
    a2(Nstrp) = fitresult.a;
    medDist2(Nstrp) = median(B2);
    
end


% pdd1 = fitdist(d1,'Normal')
% pdd2 = fitdist(d2,'Normal')
% pda = fitdist(a,'Normal')
% pdmedDist = fitdist(medDist,'Normal')
% d1ci = paramci(fitdist(d1,'Normal'))
% d2ci = paramci(fitdist(d2,'Normal'))
% aci = paramci(fitdist(a,'Normal'))


d11result = CI95(d11)
d12result = CI95(d12)
a1result = CI95(a1)
med1result = CI95(medDist1)

d21result = CI95(d21)
d22result = CI95(d22)
a2result = CI95(a2)
med2result = CI95(medDist2)

a1a2 = (a1-a2);
a1a2result = CI95(a1a2)


save('compbootstrapWithReplacement10k.mat','d11result','d12result','a1result','med1result','d21result','d22result','a2result','med2result','a1a2result')

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