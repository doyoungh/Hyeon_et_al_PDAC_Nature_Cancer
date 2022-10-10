%% Integrated clustering
% 1) Summarize membership matrix
patient_cluster=zeros(150,13);
for i=1:3
    patient_cluster(ismember(clus(:,1),strcat('C',num2str(i))),i)=1;
end
for i=1:5
    patient_cluster(ismember(clus(:,2),strcat('C',num2str(i))),i+3)=1;
end
for i=1:5
    patient_cluster(ismember(clus(:,3),strcat('C',num2str(i))),i+8)=1;
end

% 2) Perform 10,000 k-means clustering
clus_memb=[];
clus_num=[];
for i=1:10000
    if mod(i,1000)==0
        i
    end
    rng(i)
    [kmeans_idx,~]=kmeans(patient_cluster, 6);
    kmeans_idx2=zeros(150,1);
    tm_clus_num=zeros(6,1);
    for j=1:6
        if sum(patient_cluster(ismember(kmeans_idx,j),1))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),5))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),11))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=1;
            tm_clus_num(1,1)=sum(ismember(kmeans_idx,j));
        elseif sum(patient_cluster(ismember(kmeans_idx,j),2))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),6))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),11))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=2;
            tm_clus_num(2,1)=sum(ismember(kmeans_idx,j));
        elseif sum(patient_cluster(ismember(kmeans_idx,j),2))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),6))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),9))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=3;
            tm_clus_num(3,1)=sum(ismember(kmeans_idx,j));
        elseif sum(patient_cluster(ismember(kmeans_idx,j),2))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),4))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),10))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=4;
            tm_clus_num(4,1)=sum(ismember(kmeans_idx,j));
        elseif sum(patient_cluster(ismember(kmeans_idx,j),3))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),7))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),12))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=5;
            tm_clus_num(5,1)=sum(ismember(kmeans_idx,j));
        elseif sum(patient_cluster(ismember(kmeans_idx,j),3))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),8))/sum(ismember(kmeans_idx,j))>0.6 & sum(patient_cluster(ismember(kmeans_idx,j),13))/sum(ismember(kmeans_idx,j))>0.6
            kmeans_idx2(ismember(kmeans_idx,j),1)=6;
            tm_clus_num(6,1)=sum(ismember(kmeans_idx,j));
        end
    end
    clus_memb=[clus_memb;[i,kmeans_idx2']];
    clus_num=[clus_num;[i,tm_clus_num']];
end

% 3) Assign integrated clusters
clus_num_by_pat=zeros(150,7);
for i=1:150
    [clus_num_by_pat(i,:),~]=hist(clus_memb(:,i+1),[0:6]);
end

[~,int_clus]=max(clus_num_by_pat(:,2:7),[],2);

