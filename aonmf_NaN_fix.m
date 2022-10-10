function [A,S,err] = aonmf_NaN_fix(X,k,maxiter,InitA,InitS,fixedA)

[m,n] = size(X);

if isempty(InitA)
    A=zeros(m,k)+eps;
    tm=rand(sum(sum(fixedA,2)<1),k);
%     tm=tm./repmat(sum(tm,2),1,k);
    A(sum(fixedA,2)>1,1:size(fixedA,2))=fixedA(sum(fixedA,2)>1,:);
    A(sum(fixedA,2)<1,:)=tm;
    A=A./repmat(sum(A,2),1,k);
else
    A = InitA;
end

if isempty(InitS)
    S=rand(k,n);
else
    S = InitS;
end
A=normc(A);

nan_ind=isnan(X);
X(nan_ind)=0;
multiplier_1=m./sum(~nan_ind, 1);
multiplier_2=n./sum(~nan_ind, 2);

for its=1:maxiter
    S=S.*sqrt((A'*X.*multiplier_1+eps)./(A'*A*S+eps));
    A=A.*sqrt((X*S'.*multiplier_2+eps)./((A*S)*(X'*A.*multiplier_1')+eps));
%     A(sum(fixedA,2)==1,:)=[fixedA(sum(fixedA,2)==1,:),zeros(sum(sum(fixedA,2)==1),k-size(fixedA,2))];
    if its~=maxiter
%         A=A./repmat(sum(A,1),m,1);
        A=normc(A);
        err(its,1) = norm(X-A*S,'fro');
    end
    if its>1 && its~=maxiter
        if abs(err(its-1,1)-err(its,1))<eps
            break
        end
    end
end

end