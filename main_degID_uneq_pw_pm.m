function [pw,pm,OVP]=main_degID_uneq_pw_pm(trnd,wrnd,mrnd,x,y,opt)

% Compute observed statistics (t-statistic, ranksum and fold (median difference))
[tobs2,wobs2,mobs2] = stat2_uneq_ranksum(x,y);

% Compute p-values for the three statistics
pw = pval2tail(wrnd(:),wobs2); % ranksum test
pm = pval2tail(mrnd(:),mobs2); % ratio test

% Integrate the three p-values using Liptak-Stouffer's Z method
pall = [pw pm];
if opt==1
    pall = 1-normcdf(auto(-norminv(pall))); % p-value scaling and normalization
end

OVP = nwpv2(pall,3);
    
end