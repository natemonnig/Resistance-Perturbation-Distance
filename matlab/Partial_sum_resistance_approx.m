% Script for empirical testing the convergence of the partial 
% sum approximations in Theorem 10.  Used to generate Figure 4 in
% Monnig & Meyer (2016).  http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;
% explore partial sum approximations of effective resistances
N=2000;

% choose type of random graph
graphtype=3;
% type 0 = Erdos Renyi
% type 1 = two community block model
% type 2 = latent space model on unit circle
% type 3 = Watts and Strogatz small world
if graphtype==0
    % build Erdos Renyi random graph
    p=0.1;
    L=ErdosRenyi(N,p);
elseif graphtype==1
    % build stochastic block model
    split=0.5;
    Pin=0.9;
    Pout=0.005;
    
    % community membership
    membership=zeros(N,1);
    membership(rand(N,1)>split)=1;
    L=TwoCommunity(membership,Pin,Pout);
elseif graphtype==2
    % choose kernel width parameter for latent space model
    epsilon=10;

    % generate random latent locations
    theta=2*pi*rand(N,1);
    % latent space coordinates
    X=[cos(theta),sin(theta)];
    L=LatentSpace(X,epsilon);
elseif graphtype==3
    % choose number of neighbors in regular ring lattice
    % k should be even and k >> log(N)
    p=40;
    % choose rewiring probability beta
    beta=0.01;
    
    % build graph
    L=SmallWorld(N,p,beta);
end

% compute resistance matrix
Ldag=pinv(full(L));
R=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;


% pick a random pair of vertices
i=ceil(N*rand);
j=ceil(N*rand);

upperbound=zeros(N-1,2);
lowerbound=zeros(N-1,2);
% midpoint=zeros(N-1,1);
[phi,lam]=eig(full(L));
lam=diag(lam);
[lam,ind]=sort(lam);
phi=phi(:,ind);

for p=2:N
    lowerbound(p-1,1)=2/lam(N)+((phi(i,2:p)-phi(j,2:p)).^2)*(lam(2:p).^(-1)-1/lam(N));
    lowerbound(p-1,2)=2/lam(N)^2+((phi(i,2:p)-phi(j,2:p)).^2)*(lam(2:p).^(-2)-lam(N)^(-2));
    
    upperbound(p-1,1)=2/lam(p)+((phi(i,2:p)-phi(j,2:p)).^2)*(lam(2:p).^(-1)-1/lam(p));
    upperbound(p-1,2)=2/lam(p)^2+((phi(i,2:p)-phi(j,2:p)).^2)*(lam(2:p).^(-2)-lam(p)^(-2));
%     lamavg=(lam(p)+lam(N))/2;
%     midpoint(p-1)=2/lamavg+((phi(i,2:p)-phi(j,2:p)).^2)*(lam(2:p).^(-1)-1/lamavg);
end

%%

figure
subplot(1,2,1)
plot(2:N,lowerbound(:,1))
hold on
plot(2:N,upperbound(:,1))
plot(2:N,(upperbound(:,1)+lowerbound(:,1))/2)
% plot(2:N,midpoint)
plot([0 N],[lowerbound(end,1) lowerbound(end,1)])
legend('Lowerbound partial sum','Upperbound partial sum','Midpoint partial sum','\Sigma_{k=2}^n (\phi_{ki}-\phi_{kj})^2/\lambda_k')
xlabel('p')
ylabel('Approximation')
ylim([0 2*lowerbound(end,1)])

subplot(1,2,2)
plot(2:N,lowerbound(:,2))
hold on
plot(2:N,upperbound(:,2))
plot(2:N,(upperbound(:,2)+lowerbound(:,2))/2)
% plot(2:N,midpoint)
plot([0 N],[lowerbound(end,2) lowerbound(end,2)])
legend('Lowerbound partial sum','Upperbound partial sum','Midpoint partial sum','\Sigma_{k=2}^n (\phi_{ki}-\phi_{kj})^2/\lambda_k^2')
xlabel('p')
ylabel('Approximation')
ylim([0 2*lowerbound(end,2)])

figure
subplot(1,2,1)
hist(lam)
xlabel('\lambda')
ylabel('Frequency')

subplot(1,2,2)
hist(lam.^2)
xlabel('\lambda^2')
ylabel('Frequency')



figure
subplot(1,2,1)
plot(lam,1/length(lam):1/length(lam):1)
xlabel('\lambda')
ylabel('Empirical CDF')

subplot(1,2,2)
plot(lam.^2,1/length(lam):1/length(lam):1)
xlabel('\lambda^2')
ylabel('Empirical CDF')
