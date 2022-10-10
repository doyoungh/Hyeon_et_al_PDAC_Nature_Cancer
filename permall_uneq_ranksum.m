function [t,w,m] = permall_uneq_ranksum(all,nx,ny)

for i = 1:1000
    i
    ii=rand(size(all,2),1);
    [ii,jj]=sort(ii);
    xp = all(:,jj(1:nx));
    yp = all(:,jj(nx+1:nx+ny));
    [t(:,i),w(:,i),m(:,i)] = stat2_uneq_ranksum_rev2(xp,yp);
end
