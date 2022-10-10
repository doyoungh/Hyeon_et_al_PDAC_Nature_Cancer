%% Identify pairs of mRNA and protein with significant correlation
% 1) Filter overlapped genes
genes_overlap=intersect(proteins,genes);

% 2) Calculate mRNA-protein correlations
corr_protein_gene=zeros(6959,1);
for i=1:6959
    tm_gene=genes_overlap(i,1);
    tm_fc_protein=fc_merged_qnorm_prot(ismember(proteins,tm_gene),:);
    tm_fc_gene=fc_qnorm_gene(ismember(genes,tm_gene),:);
    tm_corr=zeros(size(tm_fc_protein,1),1);
    for j=1:size(tm_fc_protein,1)
        tm_corr(j,1)=corr(tm_fc_protein(j,:)',tm_fc_gene','Type','Spearman','Rows','complete');
    end
    corr_protein_gene(i,1)=max(tm_corr);
end

% 3) Generate null distributions for correlations
null_corr_protein_gene=zeros(100000,1);
for i=1:100000
    if mod(i,10000)==0;
        i
    end
    tm_ind1=randperm(6959,1);tm_entrez1=genes_overlap(tm_ind1,1);
    tm_ind2=randperm(6959,1);tm_entrez2=genes_overlap(tm_ind2,1);
    tm_fc_protein=fc_merged_qnorm_prot(ismember(proteins,tm_entrez1),:);
    tm_fc_gene=fc_qnorm_gene(ismember(genes,tm_entrez2),:);
    tm_corr=zeros(size(tm_fc_protein,1),1);
    for j=1:size(tm_fc_protein,1)
        tm_corr(j,1)=corr(tm_fc_protein(j,:)',tm_fc_gene','Type','Spearman','Rows','complete');
    end
    null_corr_protein_gene(i,1)=max(tm_corr);
end

% 4) Identify pairs of mRNA and protein with significant correlation
pval_corr_protein_gene=zeros(6959,1);
for i=1:6959
    pval_corr_protein_gene(i,1)=sum(corr_protein_gene(i,1)<null_corr_protein_gene)/100000;
end
[~,~,~,adj_pval_corr_protein_gene]=fdr_bh(pval_corr_protein_gene,0.01);
clear corr_protein_gene null_corr_protein_gene pval_corr_protein_gene

genes_sig=genes_overlap(adj_pval_corr_protein_gene<0.01,1);
clear adj_pval_corr_protein_gene

