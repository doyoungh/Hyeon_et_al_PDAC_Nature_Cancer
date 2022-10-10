%A-orthogonal NMF based on Steifel manifold
%maxiter - maximum iteration 
%k - number of basis
function [A,S] = aonmf(X,k,maxiter,InitA,InitS)
% figure out the size of data matrix
[m,N] = size(X);


% initial conditions
if isempty(InitA)
    A=rand(m,k);
else
    A = InitA;
end
if isempty(InitS)
    S=rand(k,N);
else
    S = InitS;
end


for its=1:maxiter    
    S=S.*sqrt((A'*X+eps)./(A'*A*S+eps));          
    %S=S./ repmat(sum(S,2),1,N);

    A=A.*sqrt((X*S'+eps)./((A*S)*(X'*A)+eps));      
    A=A./ repmat(sum(A,1),m,1);
    
%     err = norm(X-A*S,'fro');
%     ortho = norm(eye(size(A,2))-A'*A,'fro');
%     fprintf('iteration: %d  error %f, orthogonality of A %f \n',its,err,ortho);


end

end