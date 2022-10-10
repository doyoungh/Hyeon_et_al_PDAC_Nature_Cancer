%% 1. Normalization of mRNA expression
% 1) Identify expressed molecules
ind_exp=sum(FPKM>1,2)>196*0.5;
IDs_exp=IDs(ind_exp,1);

% 2) Normalize expression
FPKM_exp=log2(FPKM(ind_exp,:)+1);
FPKM_exp(FPKM_exp==0)=NaN;
FPKM_exp_qnorm=quantilenorm(FPKM_exp);
FPKM_exp_qnorm(isnan(FPKM_exp_qnorm))=0;
fc=FPKM_exp_qnorm-repmat(median(FPKM_exp_qnorm,2),1,196);
fc_qnorm=quantilenorm(fc);


%% 2. Normalization of protein expression
% 1) Setup parameters
PIP=70;
NUM_SET=17;
NUM_TMT=10;
IDX_PEPTIDE=1;
IDX_PIP=3;

% 2) Import scan files
cd(PATH_scans);
dir_files=dir(pwd);
dir_names={dir_files.name}';
all_names=dir_names(~ismember(dir_names,{'.','..','.DS_Store'}));
all_data=cell(NUM_SET,1);
for i=1:NUM_SET
    fid=fopen(all_names{i,1});
    temp_data=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f','delimiter','\t','HeaderLines',1);
    temp_data2=[];
    for j=1:10
        temp_data2=[temp_data2,temp_data{1,j+3}];
    end
    all_data{i,1}=table(temp_data{1,IDX_PEPTIDE}, temp_data{1,IDX_PIP}, temp_data2, 'VariableNames', {'Peptide','PIP','Intensity'});
end
clear dir_files dir_names temp_data i

% 3) Filter by PIP and log2-transform
all_data_log2=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_data=all_data{i,1};
    temp_data2=temp_data(temp_data.PIP >= PIP,:);
    temp_data3=temp_data2.Intensity;
    temp_data3(temp_data3 == 0)=NaN;
    temp_data2.Intensity=log2(temp_data3);
    all_data_log2{i,1}=temp_data2;
end
clear i temp_data temp_data2 temp_data3 all_data

% 4) Summarize scans to peptides
summary=cell(NUM_SET,1);
for i=1:NUM_SET
    disp(i/NUM_SET*100)
    temp_peptides=all_data_log2{i,1}.Peptide;
    [temp_peptides_unique,~,idx]=unique(temp_peptides);
    temp_summary_table=table(temp_peptides_unique,'VariableNames',{'Peptide'});
    for j=1:size(temp_peptides_unique,1)
        idxs=find(idx==j);
        temp_summary_table.Intensity(j,1:10)=nanmedian(all_data_log2{i,1}(idxs,:).Intensity,1);
    end
    summary{i,1}=temp_summary_table;
end
clear i j temp_peptides temp_peptides_unique idx temp_summary_table all_data_log2

% 5) Fold change calculation and within set normalization
intensity_qnorm=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_intensity=summary{i,1}.Intensity;
    intensity_qnorm{i,1}=quantilenorm(temp_intensity); 
end
fc_qnorm=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_fc=intensity_qnorm{i,1}(:,1:NUM_TMT-1)-repmat(intensity_qnorm{i,1}(:,NUM_TMT),1,NUM_TMT-1);
    fc_qnorm{i,1}=quantilenorm(temp_fc);
end
clear temp_intensity temp_fc i

% 6) Merge peptide, protein, and gene IDs
data_input=cell(NUM_SET,1);
for i=1:NUM_SET
    data_input{i,1}=table(summary{i,1}.Peptide,fc_qnorm{i,1},'VariableNames',{'Peptide','FC'});
end
clear i fc_qnorm

link_table=cell(NUM_SET,1);
data_reorder=cell(NUM_SET,1);
for i=1:NUM_SET
    link_table{i,1}=table(raw{i,1}.Peptide,pids{i,1},eids{i,1},'VariableNames',{'Peptide','Protein','EntrezID'});
    data_reorder{i,1}=outerjoin(data_input{i,1},link_table{i,1},'Type','Left','MergeKeys',true);
    data_reorder{i,1}.EntrezNum=NaN(size(data_reorder{i,1},1),1);
end
clear raw pids eids

% 7) Identify unique peptides (peptides mapped to only one gene ID)
for i=1:NUM_SET
    disp(i/NUM_SET*100);
    for j=1:size(data_reorder{i,1},1)
        if ~isempty(data_reorder{i,1}.EntrezID{j})
            data_reorder{i,1}.EntrezNum(j)=size(data_reorder{i,1}.EntrezID{j},2);
        else
            data_reorder{i,1}.EntrezNum(j)=0;
        end
    end
end
data_filtered=cell(NUM_SET,1);
for i=1:NUM_SET
    disp(i/NUM_SET*100);
    data_filtered{i,1}=data_reorder{i,1}(data_reorder{i,1}.EntrezNum==1,:);
end
clear i

% 8) Summarize peptides to proteins 
protein_data=cell(NUM_SET,1);
for i=1:NUM_SET
    disp(i/NUM_SET*100);
    T =cell2table(cell(0,size(data_filtered{i,1},2)), 'VariableNames', {'Peptide', 'FC', 'Protein', 'EntrezID', 'EntrezNum'});
    for j=1:size(data_filtered{i,1},1)
        ps=split(data_filtered{i,1}.Protein(j),'%');
        temp_T=repmat(data_filtered{i,1}(j,:),size(ps,1),1);
        for k=1:size(ps,1)
            temp_T.Protein(k)=ps(k);
        end
        T=[T; temp_T];
    end
    protein_data{i,1}=T;
end

protein_data_median=cell(NUM_SET,1);
for i=1:NUM_SET
    disp(i/NUM_SET*100);
    unique_proteins=unique(protein_data{i,1}.Protein,'stable');
    T =cell2table(cell(0,4), 'VariableNames', {'Protein', 'EntrezID', 'FC', 'Count'});
    for j=1:size(unique_proteins,1)
        index=find(strcmp(unique_proteins{j,1},protein_data{i,1}.Protein));
        if length(index) > 1
            temp_protein=protein_data{i,1}.Protein(index(1));
            temp_entrez=protein_data{i,1}.EntrezID(index(1));
            temp_FC=nanmedian(protein_data{i,1}.FC(index,:));
            temp_count=length(index);
            temp_T=table(temp_protein,temp_entrez,temp_FC,temp_count,'VariableNames', {'Protein','EntrezID','FC','Count'});
            T=[T; temp_T];
        else
            continue
        end
    end
    protein_data_median{i,1}=T;
end

% 9) Merge TMT sets
all_proteins=[];
for i=1:NUM_SET
    all_proteins=[all_proteins; protein_data_median{i,1}.Protein];
end
unique_proteins=unique(all_proteins, 'stable');
clear all_proteins

FC_merged=nan(size(unique_proteins,1), NUM_SET * (NUM_TMT-1));
for i=1:NUM_SET
    [~,b]=ismember(protein_data_median{i,1}.Protein,unique_proteins);
    FC_merged(b,(NUM_TMT-1)*(i-1)+1:(NUM_TMT-1)*i)=protein_data_median{i,1}.FC;
end

final_data=cell2table(cell(size(unique_proteins,1),2), 'VariableNames', {'Peptide', 'FC'});
final_data.Peptide=unique_proteins;
final_data.FC=FC_merged(:,ind_used==1);

% 10) Identify expressed molecules
ind_exp=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*1;  % Used for patient clustering
ind_exp2=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*0.5;  % Used for signature identification
ind_20pct=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*0.2;  % Used for mutation-protein/phosphopeptide association

% 11) Normalize expression
final_data_exp=final_data(ind_exp,:);
fc_merged_exp=final_data_exp.FC;
fc_merged_qnorm=quantilenorm2(fc_merged_exp);

final_data_exp2=final_data(ind_exp2,:);
fc_merged_exp2=final_data_exp2.FC;
proteins_union_exp2=final_data_exp2.Peptide;
fc_merged_qnorm2=quantilenorm2(fc_merged_exp2);

final_data_20pct=final_data(ind_20pct,:);
fc_merged_20pct=final_data_20pct.FC;
proteins_union_20pct=final_data_20pct.Peptide;
fc_merged_qnorm_20pct=quantilenorm2(fc_merged_20pct);


%% 3. Normalization of phosphopeptide expression
% 1) Setup parameters
PIP=70;
NUM_SET=17;
NUM_TMT=10;
IDX_PEPTIDE=1;
IDX_PIP=3;

% 2) Import scan files
cd(QUANT_PATH);
dir_files=dir(pwd);
dir_names={dir_files.name}';
all_names=dir_names(~ismember(dir_names,{'.','..','.DS_Store'}));
all_data=cell(NUM_SET,1);
for i=1:NUM_SET
    fid=fopen(all_names{i,1});
    temp_data=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f','delimiter','\t','HeaderLines',1);
    temp_data2=[];
    for j=1:10
        temp_data2=[temp_data2,temp_data{1,j+3}];
    end
    all_data{i,1}=table(temp_data{1,IDX_PEPTIDE}, temp_data{1,IDX_PIP}, temp_data2, 'VariableNames', {'Peptide','PIP','Intensity'});
end
clear dir_files dir_names temp_data i

% 3) Filter by PIP and log2-transform
all_data_log2=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_data=all_data{i,1};
    temp_data2=temp_data(temp_data.PIP >= PIP,:);
    temp_data3=temp_data2.Intensity;
    temp_data3(temp_data3 == 0)=NaN;
    temp_data2.Intensity=log2(temp_data3);
    all_data_log2{i,1}=temp_data2;
end
clear i temp_data temp_data2 temp_data3 all_data

% 4) Summarize scans to peptides
summary=cell(NUM_SET,1);
count_scans=cell(NUM_SET,1);
for i=1:NUM_SET
    disp(i/NUM_SET*100)
    temp_peptides=all_data_log2{i,1}.Peptide;
    [temp_peptides_unique,~,idx]=unique(temp_peptides);
    temp_summary_table=table(temp_peptides_unique,'VariableNames',{'Peptide'});
    for j=1:size(temp_peptides_unique,1)
        idxs=find(idx==j);
        count_scans{i,1}=[count_scans{i,1};length(idxs)];
        temp_summary_table.Intensity(j,1:10)=nanmedian(all_data_log2{i,1}(idxs,:).Intensity,1);
    end
    summary{i,1}=temp_summary_table;
end
clear i j temp_peptides temp_peptides_unique idx temp_summary_table

% 5) Fold change calculation and within set normalization
intensity_qnorm=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_intensity=summary{i,1}.Intensity;
    intensity_qnorm{i,1}=quantilenorm(temp_intensity); 
end
fc_qnorm=cell(NUM_SET,1);
for i=1:NUM_SET
    temp_fc=intensity_qnorm{i,1}(:,1:NUM_TMT-1)-repmat(intensity_qnorm{i,1}(:,NUM_TMT),1,NUM_TMT-1);
    fc_qnorm{i,1}=quantilenorm(temp_fc);
end
clear temp_intensity temp_fc i

% 6) Stack normalization
max_size=max(cellfun(@(x)size(x,1),fc_qnorm));
temp_fc=NaN(max_size*(NUM_TMT-1),NUM_SET);
for i=1:NUM_SET
    temp_fc(1:size(fc_qnorm{i,1},1)*(NUM_TMT-1),i)=fc_qnorm{i,1}(:);
end

temp_fc2=quantilenorm(temp_fc);

fc_qnorm_stacked=cell(NUM_SET,1);
for i=1:NUM_SET
    fc_qnorm_stacked{i,1}=reshape(temp_fc2(1:size(fc_qnorm{i,1},1)*(NUM_TMT-1),i),size(fc_qnorm{i,1},1),NUM_TMT-1);
end
clear temp_fc i

% 7) Merge TMT sets
temp_peptides=[];
for i=1:NUM_SET
    temp_peptides=[temp_peptides;summary{i,1}.Peptide];
end
peptides_union=unique(temp_peptides);

fc_merged=NaN(size(peptides_union,1),NUM_SET*(NUM_TMT-1));
for i=1:NUM_SET
    [~,temp_idx]=ismember(summary{i,1}.Peptide,peptides_union);
    fc_merged(temp_idx,(NUM_TMT-1)*(i-1)+1:(NUM_TMT-1)*i)=fc_qnorm_stacked{i,1};
end
clear temp_peptides temp_idx i

final_data=table(peptides_union,fc_merged(:,ind_used==1),'VariableNames', {'Peptide','FC'});

% 8) Identify expressed molecules
ind_exp=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*1;  % Used for patient clustering
ind_exp2=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*0.5;  % Used for signature identification
ind_20pct=sum(~isnan(final_data.FC),2)>=size(final_data.FC,2)*0.2;  % Used for mutation-protein/phosphopeptide association

% 9) Normalize expression
final_data_exp=final_data(ind_exp,:);
fc_merged_exp=final_data_exp.FC;
fc_merged_qnorm=quantilenorm2(fc_merged_exp);

final_data_exp2=final_data(ind_exp2,:);
fc_merged_exp2=final_data_exp2.FC;
peptides_union_exp2=final_data_exp2.Peptide;
fc_merged_qnorm2=quantilenorm2(fc_merged_exp2);

final_data_20pct=final_data(ind_20pct,:);
fc_merged_20pct=final_data_20pct.FC;
peptides_union_20pct=final_data_20pct.Peptide;
fc_merged_qnorm_20pct=quantilenorm2(fc_merged_20pct);


