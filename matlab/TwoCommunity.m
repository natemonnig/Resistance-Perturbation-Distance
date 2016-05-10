function [L] = TwoCommunity(membership,Pin,Pout)
% build graph laplacian for two community model
N=length(membership);
P=zeros(N,N);
P(membership==1,membership==1)=Pin;
P(membership==0,membership==0)=Pin;
P(membership==1,membership==0)=Pout;
P(membership==0,membership==1)=Pout;
A=zeros(N,N);
A(rand(N,N)<P)=1;
A=triu(A,1)+triu(A,1)';
L=diag(sum(A))-A;
end

