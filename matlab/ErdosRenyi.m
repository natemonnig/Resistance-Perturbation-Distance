function [L]=ErdosRenyi(N,p)
% generate Laplacian matrix for an Erdos Renyi graph
L=sprandsym(N,p);
L(L~=0)=1;
% L=spdiags(sum(A))-A;
L=spdiags(sum(L,2),0,-L);
end

