function [A,S] = aonmf_NaN(X,k,maxiter,InitA,InitS)

[m,n] = size(X);

if isempty(InitA)
    A=rand(m,k);
else
    A = InitA;
end
if isempty(InitS)
    S=rand(k,n);
else
    S = InitS;
end

nan_ind=isnan(X);
X(nan_ind)=0;
multiplier_1=m./sum(~nan_ind, 1);
multiplier_2=n./sum(~nan_ind, 2);

for its=1:maxiter
    S=S.*sqrt((A'*X.*multiplier_1+eps)./(A'*A*S+eps));
    A=A.*sqrt((X*S'.*multiplier_2+eps)./((A*S)*(X'*A.*multiplier_1')+eps));
    A=A./repmat(sum(A,1),m,1);
end

end