% Script for benchmarking speed of approximate RP-2 distance algorithm
% Used to generate Figure 3 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf

% benchmark approximation of RP2 distance
clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose size of graphs to generate
% Ns=10.^(2:4);
% Ns=2.^(7:14);
Ns=2.^(7:12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

times=zeros(length(Ns),1);
ms1=zeros(length(Ns),1);
ms2=zeros(length(Ns),1);
for i=1:length(Ns)
    disp(strcat( num2str(i),' of ',num2str(length(Ns)) ))

    % build a latent space random path with powerlaw kernel edge
    % probability:  
    % P[i~j]=p/|i-j|
    p=100;  
    N=Ns(i);
    avgdeg=p+sum(p./(p:N));
    e1=zeros(ceil(1.1*avgdeg*N),3);
    e2=zeros(ceil(1.1*avgdeg*N),3);
    count1=1;
    count2=1;
    for k=1:N
        for j=i+1:N
            if rand < (p/(j-k))
                e1(count1,:)=[k j 1];
                count1=count1+1;
            end
            if rand < (p/(j-k))
                e2(count2,:)=[k j 1];
                count2=count2+1;
            end
        end
    end
    
    disp(size(e1))
    e1(e1(:,1)==0,:)=[];
    e2(e2(:,1)==0,:)=[];
    disp(size(e1))
    
    tic
    distance_approx=drp2_approx(e1,e2);
%     distance_approx=drp2_approx_v2(e1,e2);
    times(i)=toc;
    ms1(i)=size(e1,1);
    ms2(i)=size(e2,1);
    % distance_exact=drp2_exact(e1,e2);
    
end

figure
loglog((ms1+ms2)/2,times)
hold on
loglog([1e3,1e6],[1e0,1e3])
xlabel('Number of Edges')
ylabel('Computation Time')
legend('approx. RP-2 distance','linear time')
grid on