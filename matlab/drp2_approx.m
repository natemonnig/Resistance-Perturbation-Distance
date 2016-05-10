function [ distance ] = drp2_approx(e1,e2)
% function to compute the approximate RP_2 distance
% between two graphs
%
% Uses Algorithm 2 in section 7.4 of Monnig & Meyer (2016)
% http://arxiv.org/pdf/1605.01091v1.pdf
%
% Inputs:
%   e1 = m1 x 3 matrix, encoding endoints of edges (e1(:,1),e1(:,2)) and
%       corresponding edge weights (e1(:,3)) (for symmetric graph G1)
%   e2 encodes edges of second graph G2
%
% Output:
%   Approximate RP_2 distance between G1 and G2

% user-defined error tolerance (see Theorem 9 in Monnig & Meyer, 2016)
epsilon=0.1;

% find approximate Euclidean embeddings
Z1 = Eff_Res_Approx_Embed(e1(:,1:2),e1(:,3),1e-4,epsilon);
Z2 = Eff_Res_Approx_Embed(e2(:,1:2),e2(:,3),1e-4,epsilon);

% compute distance (see Theorem 8 in Monnig & Meyer, 2016)
% if Z1 & Z2 are [k x n], and one and dd are [n x 1]
n=size(Z1,2);
one = ones(n,1);
dd = (sum(Z1.^2,1)-sum(Z2.^2,1))';
distance=sqrt(...
    2*sum(dd)^2+2*n*norm(dd,2)^2-8*(one'*Z1')*(Z1*dd)+8*(one'*Z2')*(Z2*dd) ...
    +4*norm(Z1*Z1','fro')^2+4*norm(Z2*Z2','fro')^2-8*norm(Z1*Z2','fro')^2 ...
    );
end

