function w=calc_ranksum_stat(data,ind_pat1,ind_pat2,perm)

num_pat1=sum(~isnan(data(:,ind_pat1)),2);
num_pat2=sum(~isnan(data(:,ind_pat2)),2);
num_pat_comb_uniq=unique(sort([num_pat1,num_pat2],2),'rows');
sums=cell(size(num_pat_comb_uniq,1),1);
for i=1:size(num_pat_comb_uniq,1)
    tm=num_pat_comb_uniq(i,:);
    tm2=[1:sum(tm)];
    tm_ind=zeros(perm,tm(1));
    for k=1:perm
        tm_ind(k,:)=randperm(sum(tm),tm(1));
    end
    tm_sum=zeros(perm,1);
    for j=1:perm
        tm_sum(j,1)=sum(tm2(tm_ind(j,:)));
    end
    sums{i,1}=tm_sum;
end

r=tiedrank(data');
w=zeros(size(r,2),1);
for i=1:size(r,2)
    tm_sum=sums{ismember(num_pat_comb_uniq,sort([num_pat1(i),num_pat2(i)],2),'rows'),1};
    if num_pat1(i)<=num_pat2(i)
        tm_sum_real=nansum(r(ind_pat1,i));
    else
        tm_sum_real=nansum(r(ind_pat2,i));
    end
    w(i,1)=sum(tm_sum>=tm_sum_real)/length(tm_sum);
end
end