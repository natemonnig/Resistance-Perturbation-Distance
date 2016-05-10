function [ distance ] = drp1(A1,A2)
% function to compute the exact RP_1 distance between two graphs
% 
% Definition 5 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf
%
% Input dense (symmetric) adjacency matrices

n=size(A1,1);
% generate combinatorial Laplacian matrices
L1=diag(sum(A1))-A1;
L2=diag(sum(A2))-A2;
clear A1 A2
% compute pseudoinverses
L1dag=pinv(L1);
L2dag=pinv(L2);
clear L1 L2
% compute effective resistances matrices
d1=diag(L1dag);
d2=diag(L2dag);
R1=d1*ones(1,n)+ones(n,1)*d1'-2*L1dag;
R2=d2*ones(1,n)+ones(n,1)*d2'-2*L2dag;
% find the elementwise 1-norm of the difference between resistance matrices
distance=norm(R1(:)-R2(:),1);

end