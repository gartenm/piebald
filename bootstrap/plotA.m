function plotA

load('bootstrapPVMPPMArc10k.mat', 'aresult');
A = aresult.data;
%An = A/length(A);
load('goldbootstrap410ksse.mat', 'aresult');
B = aresult.data;
%Bn = B/length(B);

figure(212);
clf;
hold on;

h1 = histogram(A);

h2 = histogram(B);

h1.Normalization = 'probability';
h1.BinWidth = 0.025;
h2.Normalization = 'probability';
h2.BinWidth = 0.025;
ylabel('normalized probability');
xlabel('relative abundance of closly apposed PVM-PPM region')


load('distroVars.mat', 'gold')
B = gold;
load('distroVars.mat', 'Avector')
A = Avector;
figure(213);
clf;
hold on;

h1 = histogram(A);

h2 = histogram(B);

h1.Normalization = 'probability';
h1.BinWidth = 2;
%h1.Values = h1.Values/max(h1.Values);

h2.Normalization = 'probability';
h2.BinWidth = 2;
ylabel('normalized probability');
xlabel('PVM-PPM distance (nm)')


%histogram([A , B],30,'Normalization','probability');