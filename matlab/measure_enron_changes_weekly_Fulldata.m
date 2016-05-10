% Script for measuring RP distance between consecutive weeks
% of the Enron email communication graph.
% Used to generate Figure 7 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;

% read in Enron email data

% Data from https://www.cs.cmu.edu/~./enron/
% processed using Leto Peel's python code, available at http://gdriv.es/letopeel
fid=fopen('data/email_links.csv');
enrondata = textscan(fid, '%d%d%s', 'delimiter',' ');
fclose(fid);
% process senders, recievers, dates
emaildate=datevec(enrondata{1,3});
senders=double(enrondata{1,1})+1;
recievers=double(enrondata{1,2})+1;
% scrap messages before year 2000
senders(emaildate(:,1)<2000)=[];
recievers(emaildate(:,1)<2000)=[];
emaildate(emaildate(:,1)<2000,:)=[];

N=max([senders;recievers]);
% aggregate emails weekly
count=1;
weeks=weeknum(datenum(emaildate));
allweeks=weeks;
allweeks(emaildate(:,1)==2001)=allweeks(emaildate(:,1)==2001)+53;
allweeks(emaildate(:,1)==2002)=allweeks(emaildate(:,1)==2002)+105;
% first Monday is Jan. 3, 2000
allmondays=7*(allweeks-1)+datenum('2000-01-03');

nweeks=max(allweeks)-min(allweeks)+1;
% store weekly effective resistance and adjacency matrices
R=zeros(N,N,nweeks);
A_all=zeros(N,N,nweeks);
    for wk=min(allweeks):max(allweeks)
        % single out desired week
        picks=allweeks==wk;
        % build and store adjacency matrix
        A=sparse([senders(picks);recievers(picks)],[recievers(picks),senders(picks)],ones(2*sum(picks),1),N,N);
        A_all(:,:,count)=full(A);
        % add a small dense matrix (to force connectivity, else R -> infty)
        A=full(A)+10^(-10)*ones(N,N);
        % build pseudoinverse and store effective resistance matrix
        Ldag=pinv(diag(sum(A))-A);
        R(:,:,count)=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;

        count=count+1;
    end
%%
% compare weekly effective resistance matrices for event detection
% (also edit distance for comparison)

RP1_dist=zeros(nweeks,1);
RP1_dist_rel=zeros(nweeks,1);
RP2_dist=zeros(nweeks,1);
RP2_dist_rel=zeros(nweeks,1);
edit_dist=zeros(nweeks,1);
edit_dist_rel=zeros(nweeks,1);
for i=2:nweeks
    RP1_dist(i)=sum(sum(abs(R(:,:,i)-R(:,:,i-1))));
    RP1_dist_rel(i)=RP1_dist(i)/sum(sum(abs(R(:,:,i-1))));
    RP2_dist(i)=norm(R(:,:,i)-R(:,:,i-1),'fro');
    RP2_dist_rel(i)=RP2_dist(i)/sum(sum(abs(R(:,:,i-1))));
    edit_dist(i)=sum(sum(abs(A_all(:,:,i)-A_all(:,:,i-1))));
    if sum(sum(abs(A_all(:,:,i-1))))>0
        edit_dist_rel(i)=edit_dist(i)/sum(sum(abs(A_all(:,:,i-1))));
    else
        edit_dist_rel(i)=0;
    end
end

mondays=min(allmondays):7:max(allmondays);

figure
hold on
for i=2:length(mondays)-1
    plot([mondays(i),mondays(i)+7],2*[RP1_dist(i),RP1_dist(i)]./max(RP1_dist),'r','LineWidth',3)
    plot([mondays(i),mondays(i)+7],2*[RP2_dist(i),RP2_dist(i)]./max(RP2_dist),'k','LineWidth',1.5)
    plot([mondays(i),mondays(i)+7],2*[edit_dist(i),edit_dist(i)]./max(edit_dist(1:end-1)),'b','LineWidth',1.5)
end
ybounds=ylim;
ylim([ybounds(1),2*ybounds(2)]);
ybounds=ylim;
xbounds=xlim;
load('data/enron_timeline.mat');
dates=cell2mat(enron_timeline(:,1))+datenum(1899,12,30);
cmap=hsv(length(dates));
deltaytext=0.0445;
for i=1:size(enron_timeline,1)
    plot([dates(i),dates(i)],[0,ybounds(2)],'Color',cmap(i,:))
end
for i=1:size(enron_timeline,1)
    text(dates(i)+4,ybounds(2)*(1-(rem(i-1,10)+1)*deltaytext),enron_timeline{i,2},...
        'EdgeColor',cmap(i,:),'FontSize',8,'BackgroundColor','w');
end
xlim([min(mondays),max(mondays)])
xticklabel_rotate(mondays(1:4:end),60,cellstr(datestr(mondays(1:4:end),2)))
xlabel('Date')
ylabel('d_{rp}(G_t,G_{t-1}) (normalized by maximum)')
box on
legend('RP1 Distance','RP2 Distance','edit distance')

figure
hold on
for i=2:length(mondays)-1
    plot([mondays(i),mondays(i)+7],[RP1_dist_rel(i),RP1_dist_rel(i)]./max(RP1_dist_rel),'r','LineWidth',3)
    plot([mondays(i),mondays(i)+7],[RP2_dist_rel(i),RP2_dist_rel(i)]./max(RP2_dist_rel),'k','LineWidth',1.5)
    plot([mondays(i),mondays(i)+7],[edit_dist_rel(i),edit_dist_rel(i)]./max(edit_dist_rel(1:end-1)),'b','LineWidth',1.5)
end
ybounds=ylim;
ylim([ybounds(1),2*ybounds(2)]);
ybounds=ylim;
xbounds=xlim;
load('data/enron_timeline.mat');
dates=cell2mat(enron_timeline(:,1))+datenum(1899,12,30);
cmap=hsv(length(dates));
deltaytext=0.0445;
for i=1:size(enron_timeline,1)
    plot([dates(i),dates(i)],[0,ybounds(2)],'Color',cmap(i,:))
end
for i=1:size(enron_timeline,1)
    text(dates(i)+4,ybounds(2)*(1-(rem(i-1,10)+1)*deltaytext),enron_timeline{i,2},...
        'EdgeColor',cmap(i,:),'FontSize',8,'BackgroundColor','w');
end
xlim([min(mondays),max(mondays)])
xticklabel_rotate(mondays(1:4:end),60,cellstr(datestr(mondays(1:4:end),2)))
xlabel('Date')
ylabel('d_{rp}(G_t,G_{t-1})/Kf(G_{t-1}) (normalized by maximum)')
box on
legend('RP1 Distance','RP2 Distance','edit distance')
