% Script for empirical testing of relationship between the magnitude
% of changes to the latent generating parameters of random graphs, and
% the resistance perturbation distance.  Used to generate Figure 6 in
% Monnig & Meyer (2016).  http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;

% choose graph size
N=2000;

% choose type of random graph
graphtype=3;
% type 1 = two community block model (change Pout)
% type 2 = latent space model on unit circle (change std dev of angle
%               changes)
% type 3 = Watts and Strogatz small world (change rewiring param beta)

nreps=50;

nparams=11;
RP1_dist=zeros(nparams,nreps);
RP2_dist=zeros(nparams,nreps);
edit_dist=zeros(nparams,nreps);

for rep_num=1:nreps
disp(strcat('Rep number = ',num2str(rep_num),' of ',num2str(nreps)))
    
if graphtype==1
    % build stochastic block model
    split=0.5;
    Pin0=0.9;
    Pout0=0.005;
    
    Pouts=Pout0:Pout0/10:2*Pout0;
%     nparams=length(Pouts);
    % community membership
    membership=zeros(N,1);
    membership(rand(N,1)>split)=1;
    L0=TwoCommunity(membership,Pin0,Pout0);
elseif graphtype==2
    % choose kernel width parameter for latent space model
    epsilon=20;
    
    % standard deviation of random angle changes
    deltas=0:0.1:1;
%     nparams=length(deltas);
    % generate random latent locations
    theta0=2*pi*rand(N,1);
    % latent space coordinates
    X=[cos(theta0),sin(theta0)];
    L0=LatentSpace(X,epsilon);
elseif graphtype==3
    % choose number of neighbors in regular ring lattice
    % k should be even and k >> log(N)
    k=40;
    % choose rewiring probability beta
    beta0=0.01;
    
    % new rewiring parameters
    betas=beta0:beta0/10:2*beta0;
%     nparams=length(betas);
    % build graph
    L0=SmallWorld(N,k,beta0);
end

% compute resistance matrix
disp(0)
Ldag=pinv(L0);
R0=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;


for l=1:nparams
    disp(l)
    
    if graphtype==1
        L1=TwoCommunity(membership,Pin0,Pouts(l));
    elseif graphtype==2
        dtheta=deltas(l)*randn(N,1);
        theta1=theta0+dtheta;
        X1=[cos(theta1),sin(theta1)];
        L1=LatentSpace(X1,epsilon);
    elseif graphtype==3
        L1=SmallWorld(N,k,betas(l));
    end
        
    Ldag=pinv(L1);
    R1=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;
    % measure Resistance Perturbation Distance
    RP1_dist(l,rep_num)=norm(R0(:)-R1(:),1);
    RP2_dist(l,rep_num)=norm(R0-R1,'fro');
    edit_dist(l,rep_num)=norm(reshape((diag(diag(L0))-L0)-(diag(diag(L1))-L1),N^2,1),1);
end

RP1_dist(:,rep_num)=RP1_dist(:,rep_num)./max(RP1_dist(:,rep_num));
RP2_dist(:,rep_num)=RP2_dist(:,rep_num)./max(RP2_dist(:,rep_num));
edit_dist(:,rep_num)=edit_dist(:,rep_num)./max(edit_dist(:,rep_num));
end

figure % error bars are standard deviation
if graphtype==1
    errorbar((Pouts-Pout0)',mean(RP1_dist,2),std(RP1_dist')')
    hold on
    errorbar((Pouts-Pout0)',mean(RP2_dist,2),std(RP2_dist')')
    errorbar((Pouts-Pout0)',mean(edit_dist,2),std(edit_dist')')
    xlabel('\Delta P_{out}')
    maxxval=max(Pouts-Pout0);
    plot([0 maxxval],[0 1],'k')
elseif graphtype==2
    errorbar(deltas',mean(RP1_dist,2),std(RP1_dist')')
    hold on
    errorbar(deltas',mean(RP2_dist,2),std(RP2_dist')')
    errorbar(deltas',mean(edit_dist,2),std(edit_dist')')
    xlabel('\sigma')
    maxxval=max(deltas);
    plot([0 maxxval],[0 1],'k')
elseif graphtype==3
    errorbar((betas-beta0)',mean(RP1_dist,2),std(RP1_dist')')
    hold on
    errorbar((betas-beta0)',mean(RP2_dist,2),std(RP2_dist')')
    errorbar((betas-beta0)',mean(edit_dist,2),std(edit_dist')')
    xlabel('\Delta \beta')
    maxxval=max(betas-beta0);
    plot([0 maxxval],[0 1],'k')
end

ylabel('d_{RP}(G^{(0)},G^{(1)}) (normalized by max)')
% legend('Resistance Perturbation Distance','Edit Distance','Perfect Correlation')
legend('RP-1 Distance','RP-2 Distance','Edit Distance','Perfect Correlation')
ylim([0 1])
xlim([0 maxxval])
box on


figure % error bars are full range of data
if graphtype==1
    errorbar((Pouts-Pout0)',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2))
    hold on
    errorbar((Pouts-Pout0)',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2))
    errorbar((Pouts-Pout0)',mean(edit_dist,2),mean(edit_dist,2)-min(edit_dist,[],2),max(edit_dist,[],2)-mean(edit_dist,2))
    xlabel('\Delta P_{out}')
    maxxval=max(Pouts-Pout0);
    plot([0 maxxval],[0 1],'k')
elseif graphtype==2
    errorbar(deltas',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2))
    hold on
    errorbar(deltas',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2))
    errorbar(deltas',mean(edit_dist,2),mean(edit_dist,2)-min(edit_dist,[],2),max(edit_dist,[],2)-mean(edit_dist,2))
    xlabel('\sigma')
    maxxval=max(deltas);
    plot([0 maxxval],[0 1],'k')
elseif graphtype==3
    errorbar((betas-beta0)',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2))
    hold on
    errorbar((betas-beta0)',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2))
    errorbar((betas-beta0)',mean(edit_dist,2),mean(edit_dist,2)-min(edit_dist,[],2),max(edit_dist,[],2)-mean(edit_dist,2))
    xlabel('\Delta \beta')
    maxxval=max(betas-beta0);
    plot([0 maxxval],[0 1],'k')
end

ylabel('d_{RP}(G^{(0)},G^{(1)}) (normalized by max)')
% legend('Resistance Perturbation Distance','Edit Distance','Perfect Correlation')
legend('RP-1 Distance','RP-2 Distance','Edit Distance','Perfect Correlation')
ylim([0 1])
xlim([0 maxxval])
box on