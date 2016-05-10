% Script for empirical testing of the relative reduction in 
% the Kirchhoff index using the result of Theorem 10 to 
% approximate Theorem 2.  Used to generate Figure 5 in
% Monnig & Meyer (2016).  http://arxiv.org/pdf/1605.01091v1.pdf

% clc; clear all; close all;

% Code for testing approximation of optimal link addition

% size of network (# vertices)
N=500;

% how many random repetitions?
nrep=50;

% compare various levels of approximation (by using subset of eigenvectors)
subnum=[N-1,1:N-2]; % first do exact, then partial sum approximations
maxdRest=zeros(N,nrep);
maxdRestGreedy=zeros(N,nrep);

% choose network type:
% 1 = Erdos Renyi
% 2 = Stochastic Block Model
% 3 = Latent Space Model on Unit Circle
% 4 = Preferential Attachment (BA model)
% 5 = Small World Model (Watts and Strogatz)

for graphtype=1:5

if graphtype==1
    graphtitle='Erdos Renyi';
elseif graphtype==2
    graphtitle='Stochastic Block Model';
elseif graphtype==3
    graphtitle='Latent Space Model';
elseif graphtype==4
    graphtitle='BA Preferential Attachment Model';
elseif graphtype==5
    graphtitle='Small World Model';
end

for rep=1:nrep
    disp(strcat({'REP '},num2str(rep),{' of '},num2str(nrep)))
    A=zeros(N,N);
    
    % build random graph
    if graphtype==1
        % choose link probability
        P=0.05;
        A(rand(N,N)<P)=1;
        A=triu(A,1)+triu(A,1)';
    elseif graphtype==2
        % choose within and between community link probabilities
        Pin=0.1;
        Pout=0.01;
        % community membership
        membership=zeros(N,1);
        membership(rand(N,1)>0.5)=1;
        % build adjacency matrices randomly from latent space models
        P=zeros(N,N);
        P(membership==1,membership==1)=Pin;
        P(membership==0,membership==0)=Pin;
        P(membership==1,membership==0)=Pout;
        P(membership==0,membership==1)=Pout;
        A(rand(N,N)<P)=1;
        A=triu(A,1)+triu(A,1)';
    elseif graphtype==3
        % build latent space around unit circle
        theta=2*pi*rand(N,1);
        epsilon=N/50;
        % latent space coordinates
        X=[cos(theta),sin(theta)];
        Xdists=zeros(N,N);
        for i=1:N
            for j=i+1:N
                Xdists(i,j)=norm(X(i,:)-X(j,:));
                Xdists(j,i)=Xdists(i,j);
            end
        end
        % use exponential connection probabilities
        P=exp(-epsilon^2*Xdists.^2);
        % build adjacency matrices randomly from latent space models
        A(rand(N,N)<P)=1;
        A=triu(A,1)+triu(A,1)';
    elseif graphtype==4
        % choose number of edges to add from each new vertex
        m=2;
        m0=2;
        % start with two connected vertices
        A(1,2)=1; A(2,1)=1;
        for i=3:N
            P=sum(A(1:i-1,1:i-1));
            P=P/sum(P);
            for j=2:i-1
                P(j)=P(j)+P(j-1);
            end
            edge1=find((rand<P),1,'first');
            edge2=find((rand<P),1,'first');
            while edge1==edge2
                edge2=find((rand<P),1,'first');
            end
            A(i,edge1)=1; A(edge1,i)=1;
            A(i,edge2)=1; A(edge2,i)=1;
        end
    elseif graphtype==5
        % choose number of neighbors in regular ring lattice
        % k should be even and k >> log(N)
        k=20;
        for i=1:N
            A(i,(i+1):min((i+k/2),N))=1;
        end
        A=A+A';
        % choose rewiring probability beta
        beta=0.1;
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
    end
    
    % experiment with optimal link addition
    [E,Lam]=eig(diag(sum(A))-A);
    Lam=diag(Lam);
    if ~issorted(Lam)
        [Lam,I]=sort(Lam);
        E=E(:,I);
    end
    
    % compute resistance matrix (for future reference)
    Ldag=pinv(diag(sum(A))-A);
    R=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;
    
    % compute change in kirchoff index for an (additional) edge
    % between each pair of vertices
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % errorbound=zeros(size(subnum));
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dR=zeros(N,N,21);
    
    k=1;
    % while (k<20 || abs(mean(maxdRest(max(2,k-5):max(k-1,2)))-maxdRest(1))>10^-3) && k<51
    while k<12
        disp(strcat(num2str(k),{' of 12'}))
        %     disp(strcat({'Using '},num2str(subnum(k)),{' eigenvectors'}))
        
        %     dR=zeros(N,N);
        
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     errorbounds=zeros(N,N);
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:N
            for j = i+1:N
                %             dR(i,j,k)=2*N*sum(((E(i,2:(subnum(k)+1))-E(j,2:(subnum(k)+1)))'./Lam(2:(subnum(k)+1))).^2)...
                %                 /(1+sum((((E(i,2:(subnum(k)+1))-E(j,2:(subnum(k)+1)))').^2)./Lam(2:(subnum(k)+1))));
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %       use partial sum bounds:
                psum1midpt=1/Lam(end)+1/Lam(subnum(k)+1)+1/2*( (E(i,2:(subnum(k)+1))-E(j,2:(subnum(k)+1))).^2 )*(2./Lam(2:(subnum(k)+1))-1/Lam(end)-1/Lam(subnum(k)+1));
                psum2midpt=1/Lam(end)^2+1/Lam(subnum(k)+1)^2+1/2*( (E(i,2:(subnum(k)+1))-E(j,2:(subnum(k)+1))).^2 )*(2./Lam(2:(subnum(k)+1)).^2-1/Lam(end)^2-1/Lam(subnum(k)+1)^2);
                dR(i,j,k)=2*N*psum2midpt/(1+psum1midpt);
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        % find indices for max estimated change
        [m,i1]=max(dR(:,:,k));
        [~,i2]=max(m);
        i1=i1(i2);
        
        % find true change by adding that edge
        A2=A;
        A2(i1,i2)=1; A2(i2,i1)=1;
        
        Ldag2=pinv(diag(sum(A2))-A2);
        R2=diag(Ldag2)*ones(1,N)+ones(N,1)*diag(Ldag2)'-2*Ldag2;
        
        maxdRest(k,rep)=norm(R(:)-R2(:),1);
        
        
        % try alternative greedy approximation method for finding optimal edge
        % to add
        dR(:,:,k)=dR(:,:,k)+dR(:,:,k)';
        i1greedy=ceil(rand*N);
        i2greedy=ceil(rand*N);
        oldi1greedy=0;
        dRrow=zeros(1,N);
        count=1;
        while i2~=oldi1greedy
            %         for j=1:N
            %             dRrow(j)=2*N*sum(((E(i2greedy,2:(subnum(k)+1))-E(j,2:(subnum(k)+1)))'./Lam(2:(subnum(k)+1))).^2)...
            %                 /(1+sum((((E(i2greedy,2:(subnum(k)+1))-E(j,2:(subnum(k)+1)))').^2)./Lam(2:(subnum(k)+1))));
            %         end
            dRrow=dR(i2greedy,:,k);
            oldi1greedy=i1greedy;
            i1greedy=i2greedy;
            [~,i2greedy]=max(dRrow);
            disp(count)
            count=count+1;
            if count==20 % re-initialize randomly
                i1greedy=ceil(rand*N);
                i2greedy=ceil(rand*N);
                oldi1greedy=0;
                count=1;
            end
            
        end
        
        % find true change by adding that edge
        A2=A;
        A2(i1greedy,i2greedy)=1; A2(i2greedy,i1greedy)=1;
        
        Ldag2=pinv(diag(sum(A2))-A2);
        R2=diag(Ldag2)*ones(1,N)+ones(N,1)*diag(Ldag2)'-2*Ldag2;
        
        maxdRestGreedy(k,rep)=norm(R(:)-R2(:),1);
        
        if maxdRestGreedy(k,rep)==maxdRest(k,rep)
            disp('Same-sies!')
        else
            disp(':(')
        end
        
        disp(maxdRest(k,rep)/maxdRest(1,rep))
        
        k=k+1;
        %     disp('maximum change in resistances:')
        %     disp(maxdR)
        %     disp('check:')
        %     disp(deltaR)
    end
    
    
end

max_corrected = maxdRest(2:(k-1),:)./repmat(max(maxdRest,[],1),10,1);
mu = mean(max_corrected,2);
mu_plus = max(max_corrected,[],2) - mu;
mu_minus = mu - min(max_corrected,[],2);

greedy_max_corrected = maxdRestGreedy(2:(k-1),:)./repmat(max(maxdRestGreedy,[],1),10,1);
greedy_mu = mean(greedy_max_corrected,2);
greedy_mu_plus = max(greedy_max_corrected,[],2) - greedy_mu;
greedy_mu_minus = greedy_mu - min(greedy_max_corrected,[],2);

figure
hold on
errorbar(subnum(2:(k-1)),mu,mu_plus,mu_minus,'r','LineWidth',3)
errorbar(subnum(2:(k-1)),greedy_mu,greedy_mu_plus,greedy_mu_minus,'b:','LineWidth',2)
hold off
xlabel('Number of Eigenvectors used in Approximation')
ylabel('\Delta KI / \Delta KI_{optimal}')
% legend('\Delta Kf / \Delta Kf_{optimal}','\Delta Kf_{greedy} / \Delta Kf_{optimal}',' \lambda_i^{-1} / \lambda_2^{-1}')
% legend('\Delta Kf / \Delta Kf_{optimal}','\Delta Kf_{greedy} / \Delta Kf_{optimal}')
legend('Exhaustive Search','Greedy Search','Location','southeast')
ylim([0 1])
xlim([0 10])
title(graphtitle)
box on
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',14)
% set(gca,'XTick',1:5)
print(graphtitle,'-depsc')
print(graphtitle,'-dpdf')
savefig(strcat(graphtitle,'.fig'))
close()
end
return
% figure
% hold on
% for i=1:N
%     for j=i+1:N
%         plot(1:size(dR,3),reshape(dR(i,j,:),1,size(dR,3)))
%     end
% end
% xlabel('Number of Eigenvectors used in Approximation')
% ylabel('Value of Approximation for pair of vertices')

figure
for m=2:9
    subplot(2,4,m-1)
    imagesc((repmat(E(:,m).^2,1,N)-2*E(:,m)*E(:,m)'+repmat((E(:,m)').^2,N,1))/Lam(m)^2)
    title(strcat('(\phi_{',num2str(m),'k}-\phi_{',num2str(m),'l})^2/\lambda_',num2str(m),'^2'))
    colorbar
end