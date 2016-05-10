function [L]=SmallWorld(N,k,beta)
% generate Laplacian matrix for a random
% Small World (Watts and Strogatz) model
A=zeros(N,N);
for i=1:k
    A=A+diag(ones(N-i,1),i)+diag(ones(i,1),i-N);
end
A=A+A';
for i=1:N
    for j=(i+1):N
        if A(i,j)==1 && rand<beta
            Aizeros=find(A(i,:)==0);
            remove=find(Aizeros==i);
            Aizeros(remove)=[];
            rewirenum=ceil(rand*length(Aizeros));
            A(i,Aizeros(rewirenum))=1; A(Aizeros(rewirenum),i)=1;
            A(i,j)=0; A(j,i)=0;
        end
    end
end
L=diag(sum(A))-A;
end

