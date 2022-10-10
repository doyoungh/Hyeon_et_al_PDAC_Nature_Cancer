function c_counts = cluster_counter(num_cluster,average_C,cg)

tree=linkage(average_C{1,num_cluster-1},'complete');
tree_count=cluster(tree,'maxclust',num_cluster);
checker=tree_count(cellfun(@str2num,cg.RowLabels));

n=hist(checker,1:num_cluster);
n2=unique(checker,'stable');
c_counts = n(n2);

end