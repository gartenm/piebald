function r = pearsonCorr(XL,YL)
% N = size(X,1)*size(X,2)*size(X,3);
% %Yn = size(Y,1)*size(Y,2);
% 
% XL = reshape(X,[N,1]);
% YL = reshape(Y,[N,1]);
N = length(XL);
mX = mean(XL);
mY = mean(YL);
%sX = sum(XL);
%xY = sum(YL);


r = (sum(XL.*YL)-N*mX*mY)/(sqrt(sum(XL.*XL)-N*mX*mX)*sqrt(sum(YL.*YL)-N*mY*mY));