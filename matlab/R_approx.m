function [ Z ] = R_approx(e,epsilon)
% code for generating approximate effective resistances
% using algorithm by Spielman and Srivastava (2008)

n = max(max(e(:,1:2)));
m = size(e,1);
% create graph Laplacian (L)
A = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[e(:,3);e(:,3)],n,n);
L = diag(sum(abs(A),2)) - A;
clear 'A' 
% create edge incidence (B) and square root of edge-weight matrices (sqrtW)
B = sparse([1:m 1:m],[e(:,1) e(:,2)],[ones(m,1) -1*ones(m,1)]);
sqrtW = sparse(1:m,1:m,e(:,3).^(1/2),m,m);
% scale = JL projection dimension
scale = ceil(log2(n))/epsilon;
% create random projection matrix Q
Q=false(scale,m);
for j=1:scale
    Q(j,:) = (rand(1,m)) > 0.5;
end
Q = Q - not(Q);
Q = Q./sqrt(scale);
% build systems of equations to solve
% Y = sparse(Q*sqrtW*B);
Y = Q*sqrtW*B;
clear Q sqrtW B
% solve each system using LAMG
inputType = 'laplacian';
solver = 'lamg';
Z=zeros(n,scale);
%---------------------------------------------------------------------
% Setup phase: construct a LAMG multi-level hierarchy
%---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
tStart  = tic;
setup   = lamg.setup(inputType, L);
tSetup  = toc(tStart);

%---------------------------------------------------------------------
% Solve phase: set up a random compatible RHS b and solve A*x=b. You can
% repeat this phase for multiple b's without rerunning the setup phase.
%---------------------------------------------------------------------
for j=1:scale
    rng(now);
    % Turn on debugging printouts during the run
    core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
%     fprintf('Solving A*x=b\n');
%     tStart = tic;
    [Z(:,j), ~, ~, ~] = lamg.solve(setup, Y(j,:)', 'errorReductionTol', 1e-8);
%     tSolve = toc(tStart);
    % Turn printouts off
    core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.INFO);
%     fprintf('------------------------------------------------------------------------\n');
end
Z=Z';
end

