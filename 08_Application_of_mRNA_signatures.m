%% 1. Normalize our data
FPKM_exp=log2(FPKM(ind_exp,:)+1);
FPKM_exp(FPKM_exp==0)=NaN;
[FPKM_exp_qnorm,mean_vals]=quantilenorm_mean(FPKM_exp);
FPKM_exp_qnorm(isnan(FPKM_exp_qnorm))=0;
fc=FPKM_exp_qnorm-repmat(median(FPKM_exp_qnorm,2),1,196);
[fc_qnorm,mean_vals_fc]=quantilenorm_mean(fc);
fc_qnorm_core=fc_qnorm(:,ind_core_mad20_c3_fix);


%% 2. Normalize data from other cohort
FPKM_other_our=FPKM_other(ismember(IDs_other,IDs_exp),:);
IDs_other_our=IDs_other(ismember(IDs_other,IDs_exp),1);
FPKM_other_our=log2(FPKM_other_our+1);
FPKM_other_our(FPKM_other_our==0)=NaN;
FPKM_other_our_qnorm=quantilenorm_given(FPKM_other_our,mean_vals);
FPKM_other_our_qnorm(isnan(FPKM_other_our_qnorm))=0;
FPKM_other_our_qnorm_reord=zeros(size(FPKM_exp_qnorm,1),size(FPKM_other_our_qnorm,2))*NaN;
for i=1:size(FPKM_exp_qnorm,1)
    if ismember(IDs_exp(i,1),IDs_other_our);
        FPKM_other_our_qnorm_reord(i,:)=FPKM_other_our_qnorm(ismember(IDs_other_our,IDs_exp(i,1)),:);
    end
end
fc_other=FPKM_other_our_qnorm_reord-repmat(median(FPKM_other_our_qnorm_reord,2),1,149);  % e.g., 149 = number of patients in other cohort
fc_other_qnorm=quantilenorm_given(fc_other,mean_vals_fc);


%% 3. Apply signatures
% 1) Calculate centroids
fc_qnorm_core_sig_other=[];
for i=1:length(gene_sig_C1_mad20_c3_fix)
    if ismember(gene_sig_C1_mad20_c3_fix(i,1),IDs_other_our);
        fc_qnorm_core_sig_other=[fc_qnorm_core_sig_other;fc_qnorm_core(ismember(IDs_exp,gene_sig_C1_mad20_c3_fix(i,1)),:)];
    end
end
for i=1:length(gene_sig_C2_mad20_c3_fix)
    if ismember(gene_sig_C2_mad20_c3_fix(i,1),IDs_other_our);
        fc_qnorm_core_sig_other=[fc_qnorm_core_sig_other;fc_qnorm_core(ismember(IDs_exp,gene_sig_C2_mad20_c3_fix(i,1)),:)];
    end
end
for i=1:length(gene_sig_C3_mad20_c3_fix)
    if ismember(gene_sig_C3_mad20_c3_fix(i,1),IDs_other_our);
        fc_qnorm_core_sig_other=[fc_qnorm_core_sig_other;fc_qnorm_core(ismember(IDs_exp,gene_sig_C3_mad20_c3_fix(i,1)),:)];
    end
end
fc_cent_other=[mean(fc_qnorm_core_sig_other(:,ismember(clus_mad20_c3_fix_core,{'C1'})),2),mean(fc_qnorm_core_sig_other(:,ismember(clus_mad20_c3_fix_core,{'C2'})),2),mean(fc_qnorm_core_sig_other(:,ismember(clus_mad20_c3_fix_core,{'C3'})),2)];

% 2) Calculate fold change for data from other cohort
fc_other_qnorm_sig=[];
for i=1:length(gene_sig_C1_mad20_c3_fix)
    if ismember(gene_sig_C1_mad20_c3_fix(i,1),IDs_other_our);
        fc_other_qnorm_sig=[fc_other_qnorm_sig;fc_other_qnorm(ismember(IDs_exp,gene_sig_C1_mad20_c3_fix(i,1)),:)];
    end
end
for i=1:length(gene_sig_C2_mad20_c3_fix)
    if ismember(gene_sig_C2_mad20_c3_fix(i,1),IDs_other_our);
        fc_other_qnorm_sig=[fc_other_qnorm_sig;fc_other_qnorm(ismember(IDs_exp,gene_sig_C2_mad20_c3_fix(i,1)),:)];
    end
end
for i=1:length(gene_sig_C3_mad20_c3_fix)
    if ismember(gene_sig_C3_mad20_c3_fix(i,1),IDs_other_our);
        fc_other_qnorm_sig=[fc_other_qnorm_sig;fc_other_qnorm(ismember(IDs_exp,gene_sig_C3_mad20_c3_fix(i,1)),:)];
    end
end

% 3) Calculate correlation with centroids
rho_other=zeros(3,149);pval_other=zeros(3,149);
[rho_other(1,:),pval_other(1,:)]=corr_pval_emp_v2(fc_cent_other(:,1),fc_qnorm_core_sig_other,fc_other_qnorm_sig,68,10000);  % e.g., 68 = number of core samples in RNA1
[rho_other(2,:),pval_other(2,:)]=corr_pval_emp_v2(fc_cent_other(:,2),fc_qnorm_core_sig_other,fc_other_qnorm_sig,58,10000);
[rho_other(3,:),pval_other(3,:)]=corr_pval_emp_v2(fc_cent_other(:,3),fc_qnorm_core_sig_other,fc_other_qnorm_sig,61,10000);
clus_other=zeros(1,149);
for i=1:149
    if sum(pval_other(:,i)<0.05)>=1;
        [~,clus_other(1,i)]=min(pval_other(:,i));
    end
end


