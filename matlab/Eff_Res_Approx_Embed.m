function [Z] = Eff_Res_Approx_Embed(e,w,tol,epsilon,type_,pfun_)
%
% calculate a Euclidean embedding for the vertices of a graph
% for which pairwise distances approximate the effective resistances 
% in the graph:
% the graph is viewed as an electrical resistive network
% where the edge weights correspond to capacitances.
% 
% Algorithm from: 
% Spielman, D. A., and Srivastava, N. Graph sparsification by effective 
% resistances. In Proceedings of the Fortieth Annual ACM Symposium on 
% Theory of Computing (2008), STOC ’08, pp. 563–568.
% 
% Details can also be found in Algorithm 1 in section 7.1 of 
% Monnig & Meyer (2016).  http://arxiv.org/pdf/1605.01091v1.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPENDENCIES:  Combinatorial Multigrid
%                available from http://www.cs.cmu.edu/~jkoutis/cmg.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   e: edges of the graph [M x 2]
%   w: weights in the graph [M x 1]
%   [tol]: specify the relative error in the linear system solution
%            1e-4 <default>
%   [epsilon]: controls accuracy in the computation of effective resistance
%              decreasing epsilon increases accuracy and computation time
%              1 <default>
%   [type_]: optional string 
%           'slm' <default>: this version use less memory.
%%%%%%%%%%%%%%%%%%ONLY USING 'spl' VERSION OF ALGORITHM%%%%%%%%%%%%%%%%%%%%
%           'spl': this implements the Spielman-Srivastava algorithm
%%%%%%%%%%%%%%%%%%ONLY USING 'spl' VERSION OF ALGORITHM%%%%%%%%%%%%%%%%%%%%
%           'org': find the exact effective resistances. Works only for
%                  small networks
%   [pfun_]: cmg_sdd(L), if already computed.
%
%
% CODE MODIFIED BY N.MONNIG TO OUTPUT FULL MATRIX Z FROM WHICH PAIRWISE
% EFFECTIVE RESISTANCES CAN BE CALCULATED
% Output:
%   Z:  embedding coordinates for vertices in space where squared euclidean
%       distances approximate effective resistances [k x n] for k=O(log n)
%       R_{ij} \approx norm(Z(:,i) - Z(:,j),2)^2
%
%
% Example usage: The path graph.
%   e = [(1:49)' (2:50)']; w = ones(length(E),1); elist = [1 50];
%   Z = Eff_Res_Approx_Embed(e,w,1e-5,1)
%
% Original code from http://www.cs.cmu.edu/~jkoutis/SpectralAlgorithms.htm
% by Richard Garcia-Lebron
%
    %% input Validation
%     [m,two] = size(e);
%     [~,two1] = size(elist);
%     if nargin > 7
    if nargin > 6
        error('Too many arguments, see help EffectiveResistancesJL.');
%     elseif nargin < 3
    elseif nargin < 2
        error('More arguments needed, see help EffectiveResistancesJL.');
%     elseif nargin == 5
    elseif nargin == 4
%         type_ = 'slm';
        type_ = 'spl';
%     elseif nargin == 4
    elseif nargin == 3
%         type_ = 'slm';
        type_ = 'spl';
        tol = 10^-4; %tolerance for the cmg solver
        epsilon = 1;
%     elseif nargin == 3
    elseif nargin == 2
%         type_ = 'slm';
        type_ = 'spl';
        tol = 10^-4; %tolerance for the cmg solver
        epsilon = 1;
%     elseif two ~= 2 ||  two1 ~= 2
%         estring = ['The first or the second input have wrong' ...
%                     ' column dimension should be [M  x 2].'];
%         error(estring);
    end

    %% Creating data for effective resitances
    numIterations = 300; %iteration for the cmg solver
    tolProb = 0.5;%use to create the johnson lindestraus projection matrix
    n = max(max(e));%number of node in the graph
    A = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[w;w],n,n); %adjacency matrix
    L = diag(sum(abs(A),2)) - A; %laplacian matrix
    clear 'A' %adjacensy matrix no needed
    m=size(e,1);
    B = sparse([1:m 1:m],[e(:,1) e(:,2)],[ones(m,1) -1*ones(m,1)]);
    [m,n] = size(B);
    W = sparse(1:length(w),1:length(w),w.^(1/2),length(w),length(w));
    scale = ceil(log2(n))/epsilon;
    %% Finding the effective resitances
    if strcmp(type_,'spl')
%         Q = sparse((rand(scale,m)) > tolProb);
        Q=false(scale,m);
        for j=1:scale
            Q(j,:) = (rand(1,m)) > tolProb;
        end
        Q = Q - not(Q);
        Q = Q./sqrt(scale);
        SYS = sparse(Q*(W)*B); % Creation the system
%         if nargin == 7 % Creating the preconditioned function
        if nargin == 6 % Creating the preconditioned function
            pfun = pfun_;
        elseif (length(L) > 600)
            pfun = cmg_sdd(L);
        end
        optimset('display','off');
        %% Solving the systems
        Z=zeros(scale,length(L));
        if (length(L) > 600) % bigger graphs
            for j=1:scale
                [Z(j,:) flag] = pcg(L,SYS(j,:)',tol,numIterations,pfun);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
            end
        else % smaller graphs
            for j=1:scale
                [Z(j,:) flag] = pcg(L,SYS(j,:)',tol,numIterations);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
            end
        end

    else
        disp('Invalid version of algorithm selected')
        Z=[];
    end
end
