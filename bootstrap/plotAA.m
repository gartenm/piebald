function plotAA

load('distroVars control.mat', 'Avector')
control = Avector;
load('distroVars KD.mat', 'Avector')
KD = Avector;


figure(312);
clf;
hold on;
h1 = histogram(control);
h2 = histogram(KD);



h1.Normalization = 'probability';
h1.BinWidth = 2;
h2.Normalization = 'probability';
h2.BinWidth = 2;
ylabel('normalized probability');
xlabel('PVM-PPM distance (nm)')

load('compbootstrapWithReplacement10k.mat', 'a1a2result')
AA = a1a2result.data;
figure(313);
clf;
hold on;

h1 = histogram(AA);
h1.Normalization = 'probability';
h1.BinWidth = 0.025;

ylabel('normalized probability');
xlabel('bootstrapped A1 - A2')
