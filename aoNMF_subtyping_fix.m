function [coph_cor,ave_C,Aa,Ss,errs]=aoNMF_subtyping_fix(data, Nbasis, Nmfiter,Numiter,clus)
% data (row: sample, column: gene)

X=data;
X=[X,-X];
X(X<0)=0;

clus_uniq=unique(clus);
fixedA=zeros(size(X,1),size(clus_uniq,2)-1)+eps;
for i=2:size(clus_uniq,2)
    fixedA(ismember(clus,clus_uniq(1,i)),i-1)=1;
end

% AONMF
Ss=cell(Numiter,1);
Aa=cell(Numiter,1);
total_C=zeros(size(data,1),size(data,1));
for k=1:Numiter
    k
    [A,S,err]=aonmf_fix(X,Nbasis,Nmfiter,[],[],fixedA);
    [~,ind]=max(A,[],2);
    for t1=1:size(A,1)
        for t2=1:size(A,1)
            if ind(t1)==ind(t2)
                C(t1,t2)=1;
            else
                C(t1,t2)=0;
            end
        end
    end                
    Ss{k}=S;
    Aa{k}=A;
    Cs{k}=C;
    total_C=[total_C+C];
    errs(:,k)=err;
end
ave_C=total_C./Numiter;
Z=linkage(ave_C,'average');
Y=pdist(ave_C);
coph_cor=cophenet(Z,Y);
end