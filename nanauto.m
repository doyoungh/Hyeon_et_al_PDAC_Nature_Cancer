function [ax,mx,stdx] = nanauto(x)
%AUTO Autoscales matrix to mean zero unit variance
%  Autoscales a matrix (x) and returns the resulting matrix (ax)
%  with mean-zero unit variance columns, a vector of means (mx) 
%  and a vector of standard deviations (stdx) used in the scaling.
%
%I/O:  [ax,mx,stdx] = auto(x);
%
%See also: MDAUTO, MDMNCN, MDRESCAL, MDSCALE, MNCN, SCALE, RESCALE

%Copyright Eigenvector Research, Inc. 1991-98
%Modified 11/93
%Checked on MATLAB 5 by BMW  1/4/97

[m,n] = size(x);
mx    = nanmean(x);
stdx  = nanstd(x);
ax    = (x-mx(ones(m,1),:))./stdx(ones(m,1),:);

