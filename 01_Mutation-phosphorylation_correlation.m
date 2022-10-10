%% Identify correlated phosphopeptides
% 1) Calculate log2 median-ratio
fc=zeros(25414,14);
for i=1:14
    fc(:,i)=nanmedian(fc_merged_qnorm(:,ind_mut(i,:)==1),2)-nanmedian(fc_merged_qnorm(:,ind_mut(i,:)==0),2);
end

% 2) Compute adjusted p-value & Determine cutoff of log2 median-ratio
pvals=zeros(25414,14*3);
mrnd_prctile=NaN(14,1);
for i=1:14
    i
    [trnd, wrnd, mrnd]=permall_uneq_ranksum(fc_merged_qnorm,sum(ind_mut(i,:)==1),sum(ind_mut(i,:)==0));
    [pw,pm,ovp]=main_degID_uneq_pw_pm(trnd,wrnd,mrnd,fc_merged_qnorm(:,ind_mut(i,:)==1),fc_merged_qnorm(:,ind_mut(i,:)==0),1);
    pvals(:,i)=pw;
    pvals(:,14+i)=pm;
    pvals(:,28+i)=ovp;
    tm=prctile(mrnd(:),[0.5,99.5]);
    mrnd_prctile(i,1)=mean(abs(tm));
end
clear trnd wrnd mrnd

% 3) Identify correlated peptides
ind_peptides_corr=logical(false(25414,14));
for i=1:14
    ind_peptides_corr(:,i)=and(fc(:,i)>mrnd_prctile(i,1),pvals(:,28+i)<0.05);
end
clear fc mrnd_prctile pvals

peptides_corr=cell(14,1);
for i=1:14
    peptides_corr{i,1}=peptides_union_exp2(ind_peptides_corr(:,i),1);
end
clear ind_peptides_corr