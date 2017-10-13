function new_wt = sirFilter(update_vec,prv_wt)
new_wt = zeros(size(prv_wt));
% clear all;close all;clc;
N=length(prv_wt);
% N = 10;
% rng('shuffle');
% pnt = rand(N,2);
% rng('shuffle');
% wt = rand(N,1);
tmp_wt = prv_wt.*update_vec;
weight = tmp_wt/sum(tmp_wt);
rng('shuffle');
delt = rand(1)/N;
cdf_val = 0;
k = 0;
for j = 0:N-1
    u = delt + j/N;
    while (u > cdf_val)
        k = k+1;
        cdf_val = cdf_val + weight(k);
    end
    rec(j+1) = k;
end
% rec;

for i = 1:N
    if histc(rec,i) ~=0
        new_wt(i) = histc(rec,i)/N;
    end
end

%     cdf_val = cdf_val + weight(i);
%     cdf_w(i) = cdf_val;
% end
% resample
% for i = 1:N
%     
% end