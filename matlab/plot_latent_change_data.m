
stdbars = figure; % error bars are standard deviation
if graphtype==1
    errorbar((Pouts-Pout0)',mean(RP1_dist,2),std(RP1_dist')','k')
    hold on
    errorbar((Pouts-Pout0)',mean(RP2_dist,2),std(RP2_dist')','r')
    errorbar((Pouts-Pout0)',mean(DeltaCon0_dist,2),std(DeltaCon0_dist')','b')
    errorbar((Pouts-Pout0)',mean(CAD_dist,2),std(CAD_dist')','g')
    xlabel('\Delta P_{out}')
    maxxval=max(Pouts-Pout0);
    plot([0 maxxval],[0 1],'y')
elseif graphtype==2
    errorbar(deltas',mean(RP1_dist,2),std(RP1_dist')','k')
    hold on
    errorbar(deltas',mean(RP2_dist,2),std(RP2_dist')','r')
    errorbar(deltas',mean(DeltaCon0_dist,2),std(DeltaCon0_dist')','b')
    errorbar(deltas',mean(CAD_dist,2),std(CAD_dist')','g')
    xlabel('\sigma')
    maxxval=max(deltas);
    plot([0 maxxval],[0 1],'y')
elseif graphtype==3
    errorbar((betas-beta0)',mean(RP1_dist,2),std(RP1_dist')','k')
    hold on
    errorbar((betas-beta0)',mean(RP2_dist,2),std(RP2_dist')','r')
    errorbar((betas-beta0)',mean(DeltaCon0_dist,2),std(DeltaCon0_dist')','b')
    errorbar((betas-beta0)',mean(CAD_dist,2),std(CAD_dist')','g')
    xlabel('\Delta \beta')
    maxxval=max(betas-beta0);
    plot([0 maxxval],[0 1],'y')
end

ylabel('d_{RP}(G^{(0)},G^{(1)}) (normalized by max)')
% legend('Resistance Perturbation Distance','Edit Distance','Perfect Correlation')
legend('RP-1 Distance','RP-2 Distance','DeltaCon_0 Distance','CAD Distance','Perfect Correlation')
ylim([0 1])
xlim([0 maxxval])
box on


maxbars = figure; % error bars are full range of data
if graphtype==1
    errorbar((Pouts-Pout0)',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2),'k')
    hold on
    errorbar((Pouts-Pout0)',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2),'r')
    errorbar((Pouts-Pout0)',mean(DeltaCon0_dist,2),mean(DeltaCon0_dist,2)-min(DeltaCon0_dist,[],2),max(DeltaCon0_dist,[],2)-mean(DeltaCon0_dist,2),'b')
    errorbar((Pouts-Pout0)',mean(CAD_dist,2),mean(CAD_dist,2)-min(CAD_dist,[],2),max(CAD_dist,[],2)-mean(CAD_dist,2),'g')
    xlabel('\Delta P_{out}')
    maxxval=max(Pouts-Pout0);
    plot([0 maxxval],[0 1],'y')
elseif graphtype==2
    errorbar(deltas',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2),'k')
    hold on
    errorbar(deltas',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2),'r')
    errorbar(deltas',mean(DeltaCon0_dist,2),mean(DeltaCon0_dist,2)-min(DeltaCon0_dist,[],2),max(DeltaCon0_dist,[],2)-mean(DeltaCon0_dist,2),'b')
    errorbar(deltas',mean(CAD_dist,2),mean(CAD_dist,2)-min(CAD_dist,[],2),max(CAD_dist,[],2)-mean(CAD_dist,2),'g')
    xlabel('\sigma')
    maxxval=max(deltas);
    plot([0 maxxval],[0 1],'y')
elseif graphtype==3
    errorbar((betas-beta0)',mean(RP1_dist,2),mean(RP1_dist,2)-min(RP1_dist,[],2),max(RP1_dist,[],2)-mean(RP1_dist,2),'k')
    hold on
    errorbar((betas-beta0)',mean(RP2_dist,2),mean(RP2_dist,2)-min(RP2_dist,[],2),max(RP2_dist,[],2)-mean(RP2_dist,2),'r')
    errorbar((betas-beta0)',mean(DeltaCon0_dist,2),mean(DeltaCon0_dist,2)-min(DeltaCon0_dist,[],2),max(DeltaCon0_dist,[],2)-mean(DeltaCon0_dist,2),'b')
    errorbar((betas-beta0)',mean(CAD_dist,2),mean(CAD_dist,2)-min(CAD_dist,[],2),max(CAD_dist,[],2)-mean(CAD_dist,2),'g')
    xlabel('\Delta \beta')
    maxxval=max(betas-beta0);
    plot([0 maxxval],[0 1],'y')
end

ylabel('d_{RP}(G^{(0)},G^{(1)}) (normalized by max)')
% legend('Resistance Perturbation Distance','Edit Distance','Perfect Correlation')
legend('RP-1 Distance','RP-2 Distance','DeltaCon_0 Distance','CAD Distance','Perfect Correlation')
ylim([0 1])
xlim([0 maxxval])
box on