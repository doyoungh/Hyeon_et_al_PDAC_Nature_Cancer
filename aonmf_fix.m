%A-orthogonal NMF based on Steifel manifold
%maxiter - maximum iteration 
%k - number of basis
function [A,S,err] = aonmf_fix(X,k,maxiter,InitA,InitS,fixedA)
% figure out the size of data matrix
[m,N] = size(X);


% initial conditions
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
% A0=A;

if isempty(InitS)
    S=rand(k,N);
else
    S = InitS;
end
% S0=S;
A=normc(A);

for its=1:maxiter
    S=S.*sqrt((A'*X+eps)./(A'*A*S+eps));          
    %S=S./ repmat(sum(S,2),1,N);
    
    A=A.*sqrt((X*S'+eps)./((A*S)*(X'*A)+eps));
%     A(sum(fixedA,2)==1,:)=[fixedA(sum(fixedA,2)==1,:),zeros(sum(sum(fixedA,2)==1),k-size(fixedA,2))];
    if its~=maxiter
%         A=A./ repmat(sum(A,1),m,1);
        A=normc(A);
        err(its,1) = norm(X-A*S,'fro');
    end
    if its>1 && its~=maxiter
        if abs(err(its-1,1)-err(its,1))<eps
            break
        end
    end
%     if mod(its,100)==0
%         ortho = norm(eye(size(A,2))-A'*A,'fro');
%         keyboard
%     end
    %     fprintf('iteration: %d  error %f, orthogonality of A %f \n',its,err,ortho);

end

end