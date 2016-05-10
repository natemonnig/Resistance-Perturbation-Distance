% Script for measuring RP distance between consecutive weeks
% of the MIT reality mining cellphone proximity graph:
% Eagle, N., and Pentland, A. Reality mining: sensing complex social 
% systems. Personal and ubiquitous computing 10, 4 (2006), 255–268.
%
% Used to generate Figure 8 in Monnig & Meyer (2016).  
% http://arxiv.org/pdf/1605.01091v1.pdf

clc; clear all; close all;

disp('Loading Reality Mining Dataset')
load('data/realitymining.mat')

% #106 has no mac ID number
s(106)=[];
% #105 is empty
s(105)=[];


%%
N=length(s);
nummeet=781032;

% first catalog all bluetooth devices of interest
blueverts=zeros(N,1);
for i=1:N
    blueverts(i)=s(i).mac;
end
% meetings = [ i , j , time(i~j) ]
meetings=zeros(nummeet,3);

count=1;
disp('Generating dynamic adjacency list from bluetooth meetings')
for i=1:N
%     disp(i)
    for j=1:length(s(i).device_macs)
        for k=1:length(s(i).device_macs{j})
            % only keep vertices of interest
            if sum(s(i).device_macs{j}(k)==blueverts)
                meetings(count,:)=[find(s(i).mac==blueverts),find(s(i).device_macs{j}(k)==blueverts),s(i).device_date(j)];
                count=count+1;
            end
        end
    end
end

clear s
clear network
%%
disp('Computing Weekly Resistance Matrices')

weeks=weeknum(meetings(:,3));
weeks(year(meetings(:,3))==2005)=weeks(year(meetings(:,3))==2005)+52;

allmondays=7*(weeks-1)+datenum('2003-12-29');

startweek=35; %first week beginning Monday Aug. 23, 2004
endweek=19+52; %last week beginning Monday May 2, 2005
numweeks=endweek-startweek+1;

R=zeros(N,N,numweeks);
A_all=zeros(N,N,numweeks);
mondays=zeros(numweeks,1);

for wk=startweek:endweek

    picks=weeks==wk;
        
    mondays(wk-startweek+1)=unique(allmondays(picks));
    
    % build weekly weighted adjacency matrix
    A=sparse(meetings(picks,1),meetings(picks,2),ones(sum(picks),1),N,N);
    A_all(:,:,wk-startweek+1)=full(A);
    A=full(A+A')+10^(-1)*ones(N,N);
        
    % compute resitance matrix
    Ldag=pinv(diag(sum(A))-A);
    R(:,:,wk-startweek+1)=diag(Ldag)*ones(1,N)+ones(N,1)*diag(Ldag)'-2*Ldag;
    
%             imagesc(R(:,:,wk))
%             drawnow
%             pause(1)
end

RP1_dist=zeros(numweeks,1);
RP1_dist_rel=zeros(numweeks,1);
RP2_dist=zeros(numweeks,1);
RP2_dist_rel=zeros(numweeks,1);
edit_dist=zeros(numweeks,1);
edit_dist_rel=zeros(numweeks,1);
for i=2:numweeks
    RP1_dist(i)=sum(sum(abs(R(:,:,i)-R(:,:,i-1))));
    RP1_dist_rel(i)=RP1_dist(i)/sum(sum(abs(R(:,:,i-1))));
    RP2_dist(i)=norm(R(:,:,i)-R(:,:,i-1),'fro');
%     RP2_dist_rel(i)=RP2_dist(i)/norm(R(:,:,i-1),'fro');
    RP2_dist_rel(i)=RP2_dist(i)/sum(sum(abs(R(:,:,i-1))));
    edit_dist(i)=sum(sum(abs(A_all(:,:,i)-A_all(:,:,i-1))));
    if sum(sum(abs(A_all(:,:,i-1))))>0
        edit_dist_rel(i)=edit_dist(i)/sum(sum(abs(A_all(:,:,i-1))));
    else
        edit_dist_rel(i)=0;
    end
end

RP1_dist_ALL=zeros(numweeks,numweeks);
for i=1:numweeks
    for j=1:numweeks
        RP1_dist_ALL(i,j)=sum(sum(abs(R(:,:,i)-R(:,:,j))));
    end
end

figure
imagesc(RP1_dist_ALL)
title('RP-1 distance between weeks, Reality Mining')
%%

figure
hold on
for i=2:length(mondays)-1
% %     plot([mondays(i),mondays(i)+7],[RP1_dist(i),RP1_dist(i)],'k','LineWidth',2)
%     plot([mondays(i),mondays(i)+7],[RP1_dist(i),RP1_dist(i)]./max(RP1_dist),'r','LineWidth',3)
%     plot([mondays(i),mondays(i)+7],[RP2_dist(i),RP2_dist(i)]./max(RP2_dist),'k','LineWidth',1.5)
    plot([mondays(i),mondays(i)+7],[edit_dist(i),edit_dist(i)]./max(edit_dist),'b','LineWidth',1.5)
end
ybounds=ylim;
ylim([0,1.5*ybounds(2)]);
ybounds=ylim;
xbounds=xlim;
load('data/reality_timeline.mat');
dates=cell2mat(reality_timeline(:,1))+datenum(1899,12,30);
cmap=hsv(length(dates));
deltaytext=0.046;
for i=1:size(reality_timeline,1)
    plot([dates(i),dates(i)],[0,ybounds(2)],'Color',cmap(i,:))
end
for i=1:size(reality_timeline,1)
    text(dates(i)+1.2,ybounds(2)*(1-(rem(i-1,5)+1)*deltaytext),reality_timeline{i,2},...
        'EdgeColor',cmap(i,:),'FontSize',8,'BackgroundColor','w');
end
xlim([min(mondays),max(mondays)])
xticklabel_rotate(mondays,60,cellstr(datestr(mondays,2)))
xlabel('Date')
% ylabel('d_{rp}(G_t,G_{t-1}) (normalized by maximum)')
ylabel('||A_t-A_{t-1}||_1 (normalized by maximum)')
box on
% legend('RP1 Distance','RP2 Distance','edit distance')
legend('edit distance')

figure
hold on
for i=2:length(mondays)-1
% %     plot([mondays(i),mondays(i)+7],[RP1_dist_rel(i),RP1_dist_rel(i)],'k','LineWidth',2)
%     plot([mondays(i),mondays(i)+7],[RP1_dist_rel(i),RP1_dist_rel(i)]./max(RP1_dist_rel),'r','LineWidth',3)
%     plot([mondays(i),mondays(i)+7],[RP2_dist_rel(i),RP2_dist_rel(i)]./max(RP2_dist_rel),'k','LineWidth',1.5)
    plot([mondays(i),mondays(i)+7],[edit_dist_rel(i),edit_dist_rel(i)]./max(edit_dist_rel),'b','LineWidth',1.5)
end
ybounds=ylim;
ylim([0,1.5*ybounds(2)]);
ybounds=ylim;
xbounds=xlim;
load('data/reality_timeline.mat');
dates=cell2mat(reality_timeline(:,1))+datenum(1899,12,30);
cmap=hsv(length(dates));
deltaytext=0.046;
for i=1:size(reality_timeline,1)
    plot([dates(i),dates(i)],[0,ybounds(2)],'Color',cmap(i,:))
end
for i=1:size(reality_timeline,1)
    text(dates(i)+1.2,ybounds(2)*(1-(rem(i-1,5)+1)*deltaytext),reality_timeline{i,2},...
        'EdgeColor',cmap(i,:),'FontSize',8,'BackgroundColor','w');
end
xlim([min(mondays),max(mondays)])
xticklabel_rotate(mondays,60,cellstr(datestr(mondays,2)))
xlabel('Date')
ylabel('d_{rp}(G_t,G_{t-1})/Kf(R_{t-1}) (normalized by maximum)')
box on
% legend('RP1 Distance','RP2 Distance','edit distance')
legend('edit distance')

% figure
% scatter(mondays(2:end),dresist(2:end))
% bounds=ylim;
% hold on
% load('data/reality_timeline.mat');
% dates=cell2mat(reality_timeline(:,1))+datenum(1899,12,30);
% cmap=hsv(length(dates));
% % cmap=cmap(randperm(size(cmap,1)),:);
% for i=1:size(reality_timeline,1)
%     plot([dates(i),dates(i)],[0,bounds(2)],'Color',cmap(i,:))
% end
% legend('d_{rp}(G_t,G_{t-1})',reality_timeline{:,2})
% % datetick('x','mm/dd/yy')
% % set(gca,'XTick',mondays)
% xlim([min(mondays),max(mondays)])
% xticklabel_rotate(mondays,60,cellstr(datestr(mondays,2)))
% xlabel('Date')
% ylabel('d_{rp}(G_t,G_{t-1})')
% 
% 
% 
% figure
% scatter(mondays(2:end),drelative(2:end))
% bounds=ylim;
% hold on
% % load('data/reality_timeline.mat');
% % dates=cell2mat(reality_timeline(:,1))+datenum(1899,12,30);
% % cmap=hsv(length(dates));
% % % cmap=cmap(randperm(size(cmap,1)),:);
% for i=1:size(reality_timeline,1)
%     plot([dates(i),dates(i)],[0,bounds(2)],'Color',cmap(i,:))
% end
% legend('d_{rp}(G_t,G_{t-1})/Kf_{t-1}',reality_timeline{:,2})
% % datetick('x','mm/dd/yy')
% % set(gca,'XTick',mondays)
% xlim([min(mondays),max(mondays)])
% xticklabel_rotate(mondays,60,cellstr(datestr(mondays,2)))
% xlabel('Date')
% ylabel('d_{rp}(G_t,G_{t-1})/Kf(G_{t-1})')
