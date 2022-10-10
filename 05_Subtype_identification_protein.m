%% 1. 1st stage ONMF clustering
% 1) Filter Top MAD genes
fc_merged_qnorm_high=fc_merged_qnorm(:,cellularity_rank<0.5451);
fc_merged_mad_high = mad(fc_merged_qnorm_high,1,2);
ind_mad_high{1,1}=find(fc_merged_mad_high>prctile(fc_merged_mad_high,90));
ind_mad_high{2,1}=find(fc_merged_mad_high>prctile(fc_merged_mad_high,80));
ind_mad_high{3,1}=find(fc_merged_mad_high>prctile(fc_merged_mad_high,70));

% 2) 1st stage ONMF clustering
for i=1:3
    i
    tm_data=fc_merged_qnorm_high(ind_mad_high{i,1},:);
    for j=2:6
        j
        [coph_cor{i}(j-1),ave_C{i}{j-1},~,~]=aoNMF_subtyping_NaN(nanauto(tm_data'),j,1000,100);
    end
end

% 3) Cluster consensus value matrix
CG_mad30_c3=clustergram(ave_C{1,3}{1,2},'cluster',3,'standardize',3,'Linkage','complete');


%% 2. 1st stage cluster assignment
% 1) Assign clusters
c_counts=cluster_counter(3,ave_C{1,3},CG_mad30_c3);
row_labels=CG_mad30_c3.RowLabels;
ind_clus_mad30_c3=cell(3,1);
for i=1:3
    ind_clus_mad30_c3{i}=cellfun(@str2num,row_labels(sum(c_counts(1:i))-c_counts(i)+1:sum(c_counts(1:i))));
end

num_patient=sum(c_counts);
clus_mad30_c3=cell(num_patient,1);
for i=1:num_patient
    for j=1:3
        if ismember(i,ind_clus_mad30_c3{j})
            clus_mad30_c3{1,i}=strcat('C',num2str(j));
        end
    end
end

% 2) Define core samples
[sil_mad30_c3,h_mad30_c3]=silhouette_NaN(nanauto(fc_merged_qnorm_high(ind_mad_high{3,1},:)'),clus_mad30_c3','distance',@sqeuclidean_NaN);


%% 3. 2nd stage ONMF clustering
% 1) Assign 1st stage ONMF clusters
clus_mad30_c3_2=cell(1,150);
clus_mad30_c3_2(1,cellularity_rank<0.5451)=clus_mad30_c3;
clus_mad30_c3_2(1,cellularity_rank>0.5451)={'C0'};

% 2) 2nd stage ONMF clustering
for i=1:3
    i
    tm_data=fc_merged_qnorm(ind_mad_high{i,1},:);
    tm_data_high=tm_data(:,cellularity_rank<0.5451);
    [~,tm_mx,tm_stdx]=auto(tm_data_high');
    tm_data_scaled=scale(tm_data',tm_mx,tm_stdx);
    for j=3:8
        j
        [coph_cor_fix{i}(j-1),ave_C_fix{i}{j-1},~,~,errs_fix{i}{j-1}]=aoNMF_subtyping_NaN_fix(tm_data_scaled,j,1000,100,clus_mad30_c3_2);
    end
end

% 3) Cluster consensus value matrix
CG_mad20_c5_fix=clustergram(ave_C_fix{1,2}{1,4},'cluster',3,'standardize',3,'Linkage','complete');


%% 4. 2nd stage cluster assignment
% 1) Extract indices for clusters
rowlabels_mad20_c5_fix=CG_mad20_c5_fix.RowLabels;
ind_clus1=cellfun(@str2num,rowlabels_mad20_c5_fix(1:33));
ind_clus2=cellfun(@str2num,rowlabels_mad20_c5_fix(34:66));
ind_clus3=cellfun(@str2num,rowlabels_mad20_c5_fix(67:88));
ind_clus4=cellfun(@str2num,rowlabels_mad20_c5_fix(89:116));
ind_clus5=cellfun(@str2num,rowlabels_mad20_c5_fix(117:end));

clus_mad20_c5_fix=cell(1,150);
for i=1:150
    if ismember(i,ind_clus1);
        clus_mad20_c5_fix{1,i}='C1';
    elseif ismember(i,ind_clus2);
        clus_mad20_c5_fix{1,i}='C2';
    elseif ismember(i,ind_clus3);
        clus_mad20_c5_fix{1,i}='C3';
    elseif ismember(i,ind_clus4);
        clus_mad20_c5_fix{1,i}='C4';
    else
        clus_mad20_c5_fix{1,i}='C5';
    end
end

% 2) Define core samples
tm_data=fc_merged_qnorm(ind_mad_high{2,1},:);
tm_data_high=tm_data(:,cellularity_rank<0.5451);
[~,tm_mx,tm_stdx]=auto(tm_data_high');
tm_data_scaled=scale(tm_data',tm_mx,tm_stdx);

[sil_mad20_c5_fix,h_mad20_c5_fix]=silhouette_NaN(tm_data_scaled,clus_mad20_c5_fix','distance',@sqeuclidean_NaN);
ind_core_mad20_c5_fix=sil_mad20_c5_fix>0;
fc_merged_qnorm2_core=fc_merged_qnorm2(:,ind_core_mad20_c5_fix);
clus_mad20_c5_fix_core=clus_mad20_c5_fix(ind_core_mad20_c5_fix);


%% 5. Identify molecular signatures
% 1) Cluster indexing
ind_Cs_mad20_c5_fix=cell(5,1);
for i=1:5
    ind_Cs_mad20_c5_fix{i}=ismember(clus_mad20_c5_fix_core,{strcat('C',num2str(i))});
end

% 2) Fold change calculation
num_peptides=size(fc_merged_qnorm2_core,1);
fc_all_mad20_c5_fix=zeros(num_peptides,5);
for i=1:5
    fc_all_mad20_c5_fix(:,i)=nanmedian(fc_merged_qnorm2_core(:,ind_Cs_mad20_c5_fix{i}),2)-nanmedian(fc_merged_qnorm2_core(:,~ind_Cs_mad20_c5_fix{i}),2);
end

% 3) empirical distribution generation
num_core_patients=size(clus_mad20_c5_fix_core,2);
core_c_counts=cellfun(@(x) sum(x),ind_Cs_mad20_c5_fix);
trnd=cell(5,1);wrnd=cell(5,1);mrnd=cell(5,1);
for i=1:5
    [trnd{i,1},wrnd{i,1},mrnd{i,1}]=permall_uneq(fc_merged_qnorm2_core,core_c_counts(i),num_core_patients-core_c_counts(i));
end

% 4) p-value calculation
pt_mad20_c5_fix=zeros(num_peptides,5);pm_mad20_c5_fix=zeros(num_peptides,5);ovp_mad20_c5_fix=zeros(num_peptides,5);
for i=1:5
    [pt_mad20_c5_fix(:,i),pm_mad20_c5_fix(:,i),ovp_mad20_c5_fix(:,i)]=main_degID_uneq(trnd{i,1},wrnd{i,1},mrnd{i,1},fc_merged_qnorm2_core(:,ind_Cs_mad20_c5_fix{i,1}),fc_merged_qnorm2_core(:,~ind_Cs_mad20_c5_fix{i,1}),1);
end
clear trnd wrnd mrnd

% 5) gene signature identification
[protein_sig_mad20_c5_fix,~]=signature_extractor(5,0,0.05,fc_all_mad20_c5_fix,ovp_mad20_c5_fix,proteins_union_exp2,fc_merged_qnorm2_core,ind_Cs_mad20_c5_fix,[]); 


