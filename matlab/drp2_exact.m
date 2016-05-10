function [ distance ] = drp2_exact(e1,e2)
% function to compute the exact RP_2 distance between two graphs
%
% Definition 5 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf
%
% Input edges:
% e1 = [m1 x 3] matrix, encoding endoints of m1 edges (e1(:,1),e1(:,2))
%       and corresponding edge weights (e1(:,3))
% e2 = [m2 x 3] ...

n=max(max(e1(:,1:2)));
% generate (dense) adjacency matrices
A1=full(sparse([e1(:,1);e1(:,2)],[e1(:,2);e1(:,1)],[e1(:,3);e1(:,3)]));
A2=full(sparse([e2(:,1);e2(:,2)],[e2(:,2);e2(:,1)],[e2(:,3);e2(:,3)]));
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
% find the Frobenius norm of the difference between resistance matrices
distance=norm(R1-R2,'fro');

end

