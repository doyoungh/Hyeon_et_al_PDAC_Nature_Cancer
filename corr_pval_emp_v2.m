function [rho,pval,null_corr]=corr_pval_emp_v2(mean_fc,fc_data1,fc_data2,num_clus_samp,perm)
% Calculate correlation between mean_fc and each column of fc_data2
% Permutation of fc_data1 and calculate mean_fc (number of 'perm' times)
% Use same null distribution for all patients in fc_data2

num_samp=size(fc_data2,2);
pval=zeros(1,num_samp);
[real_corr,~]=corr(mean_fc,fc_data2,'type','Pearson','tail','right','Rows','Pairwise');
null_corr=zeros(perm,num_samp);
for j=1:perm
    if mod(j,1000)==0;
        j
    end
    tm_rand_ind=randperm(size(fc_data1,2),num_clus_samp);
    tm_mean_fc=mean(fc_data1(:,tm_rand_ind),2);
    [null_corr(j,:),~]=corr(tm_mean_fc,fc_data2,'type','Pearson','tail','right','Rows','Pairwise');
end
rho=real_corr;

for i=1:num_samp
    pval(1,i)=sum(null_corr(:)>real_corr(1,i))/(perm*num_samp);
end

end

