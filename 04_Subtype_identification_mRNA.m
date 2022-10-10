%% 1. 1st stage ONMF clustering
% 1) Filter Top MAD genes
fc_qnorm_high=fc_qnorm(:,cellularity_rank<0.5451);  % e.g., 0.5451 = median cellularity rank
fc_qnorm_mad_high=mad(fc_qnorm_high,1,2);
ind_mad_high{1,1}=find(fc_qnorm_mad_high>prctile(fc_qnorm_mad_high,90));
ind_mad_high{2,1}=find(fc_qnorm_mad_high>prctile(fc_qnorm_mad_high,80));
ind_mad_high{3,1}=find(fc_qnorm_mad_high>prctile(fc_qnorm_mad_high,70));

% 2) 1st stage ONMF clustering
for i=1:3
    i
    tm_data=fc_qnorm_high(ind_mad_high{i,1},:);
    for j=2:6
        j
        [coph_cor{i}(j-1),ave_C{i}{j-1},~,~]=aoNMF_subtyping(auto(tm_data'),j,1000,100);
    end
end

% 3) Cluster consensus value matrix
CG_mad30_c2=clustergram(ave_C{1,3}{1,1},'cluster',3,'standardize',3,'Linkage','complete');


%% 2. 1st stage cluster assignment
% 1) Assign clusters
rowlabels_mad30_c2=CG_mad30_c2.RowLabels;
ind_clus1=cellfun(@str2num,rowlabels_mad30_c2(1:48));
ind_clus2=cellfun(@str2num,rowlabels_mad30_c2(49:end));
clus_mad30_c2=cell(1,98);
for i=1:98
    if ismember(i,ind_clus1);
        clus_mad30_c2{1,i}='C1';
    else
        clus_mad30_c2{1,i}='C2';
    end
end

% 2) Define core samples
[sil_mad30_c2,h_mad30_c2]=silhouette(auto(fc_qnorm_high(ind_mad_high{3,1},:)'),clus_mad30_c2');


%% 3. 2nd stage ONMF clustering
% 1) Assign 1st stage ONMF clusters
clus_mad30_c2_2=cell(1,196);
clus_mad30_c2_2(1,cellularity_rank<0.5451)=clus_mad30_c2;
clus_mad30_c2_2(1,cellularity_rank>0.5451)={'C0'};

% 2) 2nd stage ONMF clustering
for i=1:3
    i
    tm_data=fc_qnorm(ind_mad_high{i,1},:);
    tm_data_high=tm_data(:,cellularity_rank<0.5451);
    [~,tm_mx,tm_stdx]=auto(tm_data_high');
    tm_data_scaled=scale(tm_data',tm_mx,tm_stdx);
    for j=2:8
        j
        [coph_cor_fix{i}(j-1),ave_C_fix{i}{j-1},~,~,errs_fix{i}{j-1}]=aoNMF_subtyping_fix(tm_data_scaled,j,1000,100,clus_mad30_c2_2);
    end
end

% 3) Cluster consensus value matrix
CG_mad20_c3_fix=clustergram(ave_C_fix{1,2}{1,2},'cluster',3,'standardize',3,'Linkage','complete');


%% 4. 2nd stage cluster assignment
% 1) Extract indices for clusters
rowlabels_mad20_c3_fix=CG_mad20_c3_fix.RowLabels;
ind_clus1=cellfun(@str2num,rowlabels_mad20_c3_fix(1:70));
ind_clus2=cellfun(@str2num,rowlabels_mad20_c3_fix(71:128));
ind_clus3=cellfun(@str2num,rowlabels_mad20_c3_fix(129:end));

clus_mad20_c3_fix=cell(1,196);
for i=1:196
    if ismember(i,ind_clus1);
        clus_mad20_c3_fix{1,i}='C1';
    elseif ismember(i,ind_clus2);
        clus_mad20_c3_fix{1,i}='C2';
    else
        clus_mad20_c3_fix{1,i}='C3';
    end
end

% 2) Define core samples
tm_data=fc_qnorm(ind_mad_high{2,1},:);
tm_data_high=tm_data(:,cellularity_rank<0.5451);
[~,tm_mx,tm_stdx]=auto(tm_data_high');
tm_data_scaled=scale(tm_data',tm_mx,tm_stdx);

[sil_mad20_c3_fix,h_mad20_c3_fix]=silhouette(tm_data_scaled,clus_mad20_c3_fix');
ind_core_mad20_c3_fix=sil_mad20_c3_fix>0;
fc_qnorm_core=fc_qnorm(:,ind_core_mad20_c3_fix);
clus_mad20_c3_fix_core=clus_mad20_c3_fix(ind_core_mad20_c3_fix);


%% 5. Identify molecular signatures
% 1) Cluster indexing
ind_C1=ismember(clus_mad20_c3_fix_core,{'C1'});
ind_C2=ismember(clus_mad20_c3_fix_core,{'C2'});
ind_C3=ismember(clus_mad20_c3_fix_core,{'C3'});

% 2) Fold change calculation
fc_all_mad20_c3_fix=zeros(12967,3);
fc_all_mad20_c3_fix(:,1)=median(fc_qnorm_core(:,ind_C1),2)-median(fc_qnorm_core(:,~ind_C1),2);
fc_all_mad20_c3_fix(:,2)=median(fc_qnorm_core(:,ind_C2),2)-median(fc_qnorm_core(:,~ind_C2),2);
fc_all_mad20_c3_fix(:,3)=median(fc_qnorm_core(:,ind_C3),2)-median(fc_qnorm_core(:,~ind_C3),2);

% 3) empirical distribution generation
[trnd1,wrnd1,mrnd1]=permall_uneq(fc_qnorm_core,68,187-68);  % e.g., 187 = number of core samples, 68 = number of core samples in RNA1
[trnd2,wrnd2,mrnd2]=permall_uneq(fc_qnorm_core,58,187-58);
[trnd3,wrnd3,mrnd3]=permall_uneq(fc_qnorm_core,61,187-61);

% 4) p-value calculation
pt_mad20_c3_fix=zeros(12967,3);pm_mad20_c3_fix=zeros(12967,3);ovp_mad20_c3_fix=zeros(12967,3);
[pt_mad20_c3_fix(:,1),pm_mad20_c3_fix(:,1),ovp_mad20_c3_fix(:,1)]=main_degID_uneq(trnd1,wrnd1,mrnd1,fc_qnorm_core(:,ind_C1),fc_qnorm_core(:,~ind_C1),1);
[pt_mad20_c3_fix(:,2),pm_mad20_c3_fix(:,2),ovp_mad20_c3_fix(:,2)]=main_degID_uneq(trnd2,wrnd2,mrnd2,fc_qnorm_core(:,ind_C2),fc_qnorm_core(:,~ind_C2),1);
[pt_mad20_c3_fix(:,3),pm_mad20_c3_fix(:,3),ovp_mad20_c3_fix(:,3)]=main_degID_uneq(trnd3,wrnd3,mrnd3,fc_qnorm_core(:,ind_C3),fc_qnorm_core(:,~ind_C3),1);
clear trnd1 wrnd1 mrnd1 trnd2 wrnd2 mrnd2 trnd3 wrnd3 mrnd3

% 5) gene signature identification
gene_sig_C1_mad20_c3_fix=IDs_exp(fc_all_mad20_c3_fix(:,1)>0.58 & ovp_mad20_c3_fix(:,1)<0.06 & median(fc_qnorm_core(:,ind_C1),2)>0 & median(fc_qnorm_core(:,~ind_C1),2)<0 & median(fc_qnorm_core(:,ind_C1),2)>median(fc_qnorm_core(:,ind_C2),2) & median(fc_qnorm_core(:,ind_C1),2)>median(fc_qnorm_core(:,ind_C3),2));
gene_sig_C2_mad20_c3_fix=IDs_exp(fc_all_mad20_c3_fix(:,2)>0.58 & ovp_mad20_c3_fix(:,2)<0.06 & median(fc_qnorm_core(:,ind_C2),2)>0 & median(fc_qnorm_core(:,~ind_C2),2)<0 & median(fc_qnorm_core(:,ind_C2),2)>median(fc_qnorm_core(:,ind_C1),2) & median(fc_qnorm_core(:,ind_C2),2)>median(fc_qnorm_core(:,ind_C3),2));
gene_sig_C3_mad20_c3_fix=IDs_exp(fc_all_mad20_c3_fix(:,3)>0.58 & ovp_mad20_c3_fix(:,3)<0.06 & median(fc_qnorm_core(:,ind_C3),2)>0 & median(fc_qnorm_core(:,~ind_C3),2)<0 & median(fc_qnorm_core(:,ind_C3),2)>median(fc_qnorm_core(:,ind_C1),2) & median(fc_qnorm_core(:,ind_C3),2)>median(fc_qnorm_core(:,ind_C2),2));
gene_sig_mad20_c3_fix=[[gene_sig_C1_mad20_c3_fix,repmat({'C1'},length(gene_sig_C1_mad20_c3_fix),1)];[gene_sig_C2_mad20_c3_fix,repmat({'C2'},length(gene_sig_C2_mad20_c3_fix),1)];[gene_sig_C3_mad20_c3_fix,repmat({'C3'},length(gene_sig_C3_mad20_c3_fix),1)]];



