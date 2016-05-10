% Script for testing the scaling of the DeltaCon distance vs. the 
% resistance perturbation distance for simple graphs.
% Used to generate Figure 2 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;

Ns=2.^(4:9);

complete_RP=zeros(size(Ns));
complete_DC=zeros(size(Ns));

star_RP=zeros(size(Ns));
star_DC=zeros(size(Ns));

path_RP=zeros(size(Ns));
path_DC=zeros(size(Ns));

for i=1:length(Ns)
    
    % complete graph adjacency matrix
    A1=ones(Ns(i),Ns(i))-eye(Ns(i));
    A2=A1; A2(1,2)=2; A2(2,1)=2;
    % measure distances
    complete_RP(i)=drp1(A1,A2);
    complete_DC(i)=deltacon0(A1,A2);
    
    % star graph adjacency matrix
    A1=zeros(Ns(i),Ns(i));
    A1(2:end,1)=1; A1(1,2:end)=1;
    A2=A1; A2(1,2)=2; A2(2,1)=2;
    % measure distances
    star_RP(i)=drp1(A1,A2);
    star_DC(i)=deltacon0(A1,A2);
    
    % path graph adjacency matrix (turn into a cycle)
    A1=diag(ones(Ns(i)-1,1),1)+diag(ones(Ns(i)-1,1),-1);
    A2=A1; A2(1,end)=1; A2(end,1)=1;
    % measure distances
    path_RP(i)=drp1(A1,A2);
    path_DC(i)=deltacon0(A1,A2);
    
end


figure
subplot(1,2,1)
loglog(Ns,complete_RP,'-ko','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
hold on
loglog(Ns,star_RP,'-k*','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
loglog(Ns,path_RP,'-ks','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
legend('complete graph','star graph','path graph')
title('Resistance Perturbation Distance')
xlabel('order of graph (n=|V|)')
ylabel('d_{rp1}(G,G+e)')
grid on


subplot(1,2,2)
loglog(Ns,complete_DC,'-ko','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
hold on
loglog(Ns,star_DC,'-k*','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
loglog(Ns,path_DC,'-ks','LineWidth',0.1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','k',...
                       'MarkerSize',10)
legend('complete graph','star graph','path graph')
title('DeltaCon_0 Distance')
xlabel('order of graph (n=|V|)')
ylabel('d_{DC0}(G,G+e)')
grid on
