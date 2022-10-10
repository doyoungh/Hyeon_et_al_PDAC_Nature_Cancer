function d2=sqeuclidean_NaN(XI,XJ)
% distance function for executing pdist
% Euclidean distance accepting NaNs

num_row=size(XJ,1);
d2=zeros(num_row,1);
for i=1:num_row
    cnt=sum(~isnan(XI)&~isnan(XJ(i,:)));
    d2(i,1)=nansum((XI-XJ(i,:)).^2)*(length(XI)/cnt);
end
end