function [coph_cor,ave_C,Aa,Ss]=aoNMF_subtyping_NaN(data, Nbasis, Nmfiter,Numiter)
% data (row: sample, column: gene)

X=data;
X=[X,-X];
X(X<0)=0;

% AONMF
Ss=cell(Numiter,1);
Aa=cell(Numiter,1);
total_C=zeros(size(data,1),size(data,1));
for k=1:Numiter
    [A,S]=aonmf_NaN(X,Nbasis,Nmfiter,[],[]);
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
end
ave_C=total_C./Numiter;
Z=linkage(ave_C,'average');
Y=pdist(ave_C);
coph_cor=cophenet(Z,Y);
end