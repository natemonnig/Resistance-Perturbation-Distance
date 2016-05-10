function [L]=LatentSpace(X,epsilon)
% build random network from latent space model on unit circle
% use exponential connection probabilities
N=size(X,1);
P=exp(-epsilon^2*squareform(pdist(X).^2));
% build adjacency matrices randomly from latent space models
A=zeros(N,N);
A(2*rand(N,N)<P)=1;
A=max(A,A');
L=diag(sum(A))-A;
end

