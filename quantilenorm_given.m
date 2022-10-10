function normdata=quantilenorm_given(data,mean_vals)
% Quantile normalization using given quantile vectors

normdata=data;
dataSize=size(data);
nanvals=isnan(data);
numNans=sum(nanvals);
ndx=ones(dataSize);
N_quantile=length(mean_vals);
N=dataSize(1);
rr=cell(dataSize(2),1);

for col=1:dataSize(2)
    [~,ndx(:,col)]=sort(data(:,col));
    rr{col}=sort(tiedrank(data(~nanvals(:,col),col)));
end

for col=1:dataSize(2)
    M=N-numNans(col);
    normdata(ndx(1:M,col),col)=interp1(1:N_quantile,mean_vals,1+((N_quantile-1)*(rr{col}-1)/(M-1)));
end

end