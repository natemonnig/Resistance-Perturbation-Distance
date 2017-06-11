clear all; close all;
n = 200;
p = 0.1:0.01:0.2;
n_reps = 100;

deltacon_dist=zeros(size(p));
RP1_dist=zeros(size(p));

for i = 1:length(p)
    disp(['Running #' num2str(i) ' of ' num2str(length(p))])
    for k = 1:n_reps
        A1 = full(sprandsym(n,p(i)));
        A1(A1 ~= 0) = 1;
        A1 = A1-diag(diag(A1));
        % just remove the first edge in A, I think this is "random" enough...
        I = find(A1,1,'first');
        A2 = A1;
        A2(I) = 0;
        RP1_dist(i) = RP1_dist(i) + drp1(A1,A2);
        deltacon_dist(i) = deltacon_dist(i) + deltacon0(A1,A2);
    end
end

figure(1)
hold on
semilogy(p,RP1_dist./max(RP1_dist),'ro')
semilogy(p,deltacon_dist./max(deltacon_dist),'bs')
legend('RP1 Distance','DeltaCon_0 Distance')
xlabel('edge probability')
ylabel('mean d(G,G-e_1) (normalized by maximum)')
title('Erdos-Renyi Edge Submodularity')