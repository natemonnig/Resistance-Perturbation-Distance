function [ d_DC0, S1, S2 ] = deltacon0(A1,A2)
% function to compute the DeltaCon_0 distance between two graphs
% input dense adjacency matrices:
% eqn. 3.3 in Koutra et al. (2013) http://arxiv.org/pdf/1304.4657v1.pdf

n=size(A1,1);

D1=sum(A1);
D2=sum(A2);

epsilon1=1/(1+max(D1));
epsilon2=1/(1+max(D2));

S1=inv(eye(n)+epsilon1^2*diag(D1)-epsilon1*A1);
S2=inv(eye(n)+epsilon2^2*diag(D2)-epsilon2*A2);

d_DC0=norm(S1.^(1/2)-S2.^(1/2),'fro');

end

