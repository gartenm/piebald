function [T1,T2] = costes(C1,C2)
abFit = polyfit(C1,C2,1);
a = abFit(1);
b = abFit(2);
figure(532);
clf;
hold on;
plot(C1,C2,'.');
xlabel('Channel 1 intensity');
ylabel('Channel 2 intensity');

Xrange = [1:nanmax(C1)];
plot(Xrange,a*Xrange+b);

figure(4217);
clf;
hold on;

fun = @(T)RT(C1,C2,a,b,T);
T0 = 0.1*nanmax(C1);
T = fminsearch(fun,T0);


T1 = T;
T2 = a*T+b;


function r = RT(C1,C2,a,b,T)
T1 = T;
T2 = a*T+b;

%idx = ((C1<T1)&(C2<T2));
idx = ((C1<T1)|(C2<T2));

r = abs(pearsonCorr(C1(idx),C2(idx)));
%r = (pearsonCorr(C1(idx),C2(idx)));
figure(4217);
plot(T,r,'s')