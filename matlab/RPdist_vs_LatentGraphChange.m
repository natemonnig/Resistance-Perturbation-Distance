% Script for empirical testing of relationship between the magnitude
% of changes to the latent generating parameters of random graphs, and
% the resistance perturbation distance.  Used to generate Figure 6 in
% Monnig & Meyer (2016).  http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;

% choose graph size
N=2000;

% choose type of random graph
graphtype=1;
% type 1 = two community block model (change Pout)
% type 2 = latent space model on unit circle (change std dev of angle
%               changes)
% type 3 = Watts and Strogatz small world (change rewiring param beta)

nreps=50;

nparams=11;
RP1_dist=zeros(nparams,nreps);
RP2_dist=zeros(nparams,nreps);
DeltaCon0_dist=zeros(nparams,nreps);
CAD_dist=zeros(nparams,nreps);

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
A0=diag(diag(L0))-L0;

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
    A1=diag(diag(L1))-L1;
    % measure Resistance Perturbation Distance
    RP1_dist(l,rep_num)=norm(R0(:)-R1(:),1);
    RP2_dist(l,rep_num)=norm(R0-R1,'fro');
    DeltaCon0_dist(l,rep_num)=deltacon0(A0,A1);
    CAD_dist(l,rep_num)=CAD_distance(A0,A1,R0,R1);
end

RP1_dist(:,rep_num)=RP1_dist(:,rep_num)./max(RP1_dist(:,rep_num));
RP2_dist(:,rep_num)=RP2_dist(:,rep_num)./max(RP2_dist(:,rep_num));
DeltaCon0_dist(:,rep_num)=DeltaCon0_dist(:,rep_num)./max(DeltaCon0_dist(:,rep_num));
CAD_dist(:,rep_num)=CAD_dist(:,rep_num)./max(CAD_dist(:,rep_num));
end

clear A0 A1 L0 L1 Ldag R0 R1

% run plotting script
plot_latent_change_data

% save figures and data
% Note:  data can be reloaded, and plotted without re-running this script
if ~isdir('figures')
    mkdir('figures')
end

if graphtype==1
    saveas(stdbars,'figures/RP_dist_vs_two_community_stdbars.fig')
    saveas(maxbars,'figures/RP_dist_vs_two_community_maxbars.fig')
    save('figures/RP_dist_two_community_vars.mat')
elseif graphtype==2
    saveas(stdbars,'figures/RP_dist_vs_latent_circle_stdbars.fig')
    saveas(maxbars,'figures/RP_dist_vs_latent_circle_maxbars.fig')
    save('figures/RP_dist_latent_circle_vars.mat')
elseif graphtype==3
    saveas(stdbars,'figures/RP_dist_vs_small_world_stdbars.fig')
    saveas(maxbars,'figures/RP_dist_vs_small_world_maxbars.fig')
    save('figures/RP_dist_small_world_vars.mat')
end