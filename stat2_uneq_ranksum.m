function [t,w,mdiff] = stat2_uneq_ranksum(x,y)
m = size(x,1); nx = sum(~isnan(x),2); ny = sum(~isnan(y),2);
nx2 = size(x,2); ny2 = size(y,2);

% Median difference
mdiff = nanmedian(x,2)-nanmedian(y,2);

% two-sample t-test (unequal variance of two groups is assumed => default
% setting of ttest2.m)

difference = nanmean(x,2) - nanmean(y,2);
% dfe = nx + ny - 2;
% sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
% se = sPooled .* sqrt(1./nx + 1./ny);

s2x = nanvar(x,[],2);
s2y = nanvar(y,[],2);

s2xbar = s2x ./ nx;
s2ybar = s2y ./ ny;
dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
se = sqrt(s2xbar + s2ybar);

%keyboard
%nnk=prctile(se,25); se(find(se<nnk))=nnk; %standard deviation ceiling
t = difference ./ se;

% Ranksum statistic
w=calc_ranksum_stat([x y],[1:nx2],[nx2+1:nx2+ny2],1000);

t(find(isnan(t)))=0;
t(find(isinf(t)))=0;

end
