%% Identify associations
% 1) Calculate statistics
onesamp_p=cell(24,2);
onesamp_t=cell(24,2);
med_sig=cell(24,2);
for i=1:24
    tm_ind_mut=Mut_type_index(:,i)>0;
    tm_ind_nomut=Mut_type_index(:,i)==0;
    tm_ind_prot=ismember(genes_20pct_prot,freq_genes(i,1));
    tm_ind_phos=ismember(genes_20pct_phos,freq_genes(i,1));
    if sum(tm_ind_prot)~=0
        tm_fc_prot_mut=fc_merged_qnorm_20pct_prot(tm_ind_prot,tm_ind_mut);
        tm_fc_prot_nomut=fc_merged_qnorm_20pct_prot(tm_ind_prot,tm_ind_nomut);
        [~,p_left,~,stats]=ttest(tm_fc_prot_mut',0,'Tail','left');
        [~,p_right,~,~]=ttest(tm_fc_prot_mut',0,'Tail','right');
        onesamp_p{i,1}=[p_left;p_right]';onesamp_t{i,1}=stats.tstat';
        med_sig{i,1}=nanmedian(tm_fc_prot_mut,2)>prctile(tm_fc_prot_nomut',70)' | nanmedian(tm_fc_prot_mut,2)<prctile(tm_fc_prot_nomut',30)';
    end
    if sum(tm_ind_phos)~=0
        tm_fc_phos_mut=fc_merged_qnorm_20pct_phos(tm_ind_phos,tm_ind_mut);
        tm_fc_phos_nomut=fc_merged_qnorm_20pct_phos(tm_ind_phos,tm_ind_nomut);
        [~,p_left,~,stats]=ttest(tm_fc_phos_mut',0,'Tail','left');
        [~,p_right,~,~]=ttest(tm_fc_phos_mut',0,'Tail','right');
        onesamp_p{i,2}=[p_left;p_right]';onesamp_t{i,2}=stats.tstat';
        med_sig{i,2}=nanmedian(tm_fc_phos_mut,2)>prctile(tm_fc_phos_nomut',70)' | nanmedian(tm_fc_phos_mut,2)<prctile(tm_fc_phos_nomut',30)';
    end
end

% 2) Identify associations
ind_assoc=cell(24,2);
for i=1:24
    for j=1:2
        if ~isempty(onesamp_p{i,j})
            tm=zeros(size(onesamp_p{i,j},1),1);
            tm(onesamp_p{i,j}(:,2)<0.05 & med_sig{i,j}==1 & onesamp_t{i,j}>0)=1;
            tm(onesamp_p{i,j}(:,1)<0.05 & med_sig{i,j}==1 & onesamp_t{i,j}<0)=-1;
            ind_assoc{i,j}=tm;
        end
    end
end
clear onesamp_p onesamp_t med_sig

proteins_assoc=cell(24,2);
peptides_assoc=cell(24,2);
for i=1:24
    tm_prot=proteins_union_20pct(ismember(genes_20pct_prot,freq_genes(i,1)),1);
    tm_phos=peptides_union_20pct(ismember(genes_20pct_phos,freq_genes(i,1)),1);
    if ~isempty(ind_assoc{i,1})
        if sum(ind_assoc{i,1}~=0)>0
            proteins_assoc{i,1}=tm_prot(ind_assoc{i,1}>0,1);
            proteins_assoc{i,2}=tm_prot(ind_assoc{i,1}<0,1);
        end
    end
    if ~isempty(ind_assoc{i,2})
        if sum(ind_assoc{i,2}~=0)>0
            peptides_assoc{i,1}=tm_phos(ind_assoc{i,2}>0,1);
            peptides_assoc{i,2}=tm_phos(ind_assoc{i,2}<0,1);
        end
    end
end
clear ind_assoc



