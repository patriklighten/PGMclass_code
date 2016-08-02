% MHUNIFORMTRANS
%
%  MCMC Metropolis-Hastings transition function that
%  utilizes the uniform proposal distribution.
%  A - The current joint assignment.  This should be
%      updated to be the next assignment
%  G - The network
%  F - List of all factors
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function A = MHUniformTrans(A, G, F)

% Draw proposed new state from uniform distribution
A_prop = ceil(rand(1, length(A)) .* G.card);

p_acceptance = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
% Compute acceptance probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%M = ComputeApproxMarginalsBP(F,[]);%compute the stationary distribution 

%compute the stationary distribution for two states respectively. 
Pi_1=0;
Pi_2=0;
%for i=1:length(A)
   % Pi_1=Pi_1*M(i).val(A(i));
   % Pi_2=Pi_1*M(i).val(A_prop(i));
%end

%becuse the pi is uniform Q(x->x')==Q(x'->x)
V=1:length(A);
factorIndx=unique([G.var2factors{V}]);
%becuse the pi is uniform Q(x->x')==Q(x'->x)

%for i=1:length(factorIndx)
  %  Pi_1=Pi_1+log(GetValueOfAssignment(F(factorIndx(i)),A(:,F(factorIndx(i)).var)));
  %  Pi_2=Pi_2+log(GetValueOfAssignment(F(factorIndx(i)),A_prop(:,F(factorIndx(i)).var)));
%end
Pi_1 = LogProbOfJointAssignment(F, A);
Pi_2=LogProbOfJointAssignment(F, A_prop);

p_acceptance=min(1,exp(Pi_2-Pi_1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accept or reject proposal
if rand() < p_acceptance
    % disp('Accepted');
    A = A_prop;
end