function result = CI95 (data)
interval = 95;
cut = (100-interval)/200;
N = length(data);
cutoff = round(cut*N);
B = sort(data);
low = B(cutoff+1);
up = B(end-cutoff);
center95 = data((data>=low) & (data<=up));

% center95 = B(cutoff:(end-cutoff));
% result.low = center95(1);
% result.up = center95(end);

result.low = low;
result.up = up;


result.datamean = mean(data);
result.datamedian = median(data);
result.data = data;
result.center95 = center95;