function getAerror

load('goldbootstrap410ksse.mat', 'aresult')
Agold = aresult.data;

load('bootstrapPVMPPMArc10k.mat', 'aresult')
APPMPVM = aresult.data;


Rclose = Agold./APPMPVM;
Rlumen = (1-Agold)./(1-APPMPVM);

%RR = Rclose./Rlumen;

Rclose = CI95(Rclose)

Rlumen = CI95(Rlumen)

RR = (Rclose.center95)./(Rlumen.center95);

%RR = CI95(RR)

hist(RR)
mean(RR)
median(RR)
min(RR)
max(RR)