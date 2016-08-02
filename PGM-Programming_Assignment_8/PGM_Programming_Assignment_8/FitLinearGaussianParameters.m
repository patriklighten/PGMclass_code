function [Beta sigma] = FitLinearGaussianParameters(X, U)

% Estimate parameters of the linear Gaussian model:
% X|U ~ N(Beta(1)*U(1) + ... + Beta(n)*U(n) + Beta(n+1), sigma^2);

% Note that Matlab/Octave index from 1, we can't write Beta(0).
% So Beta(n+1) is essentially Beta(0) in the text book.

% X: (M x 1), the child variable, M examples
% U: (M x N), N parent variables, M examples
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

M = size(U,1);
N = size(U,2);

Beta = zeros(N+1,1);
sigma = 1;

% collect expectations and solve the linear system
% A = [ E[U(1)],      E[U(2)],      ... , E[U(n)],      1     ; 
%       E[U(1)*U(1)], E[U(2)*U(1)], ... , E[U(n)*U(1)], E[U(1)];
%       ...         , ...         , ... , ...         , ...   ;
%       E[U(1)*U(n)], E[U(2)*U(n)], ... , E[U(n)*U(n)], E[U(n)] ]

% construct A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[E_U,V_U]=FitGaussianParameters(U);%E_U is 1*N matrix, each of which entry is the mean for a certain variables 

A=ones(N+1,N+1);
A(1,:)=[E_U,1];

cov_Uij=zeros(N,N);
for i=2:N+1
    
    [miu,V_XpluesY]=FitGaussianParameters(U+repmat(U(:,i-1),1,N));
    
    cov_Uij(i-1,:)=(V_XpluesY.^2-V_U.^2-V_U(i-1)^2)/2;
    E_XY=cov_Uij(i-1,:)+E_U*E_U(i-1);
    
    A(i,:)=[E_XY,E_U(i-1)];
end

    


% B = [ E[X]; E[X*U(1)]; ... ; E[X*U(n)] ]
% construct B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[E_X,V_X]=FitGaussianParameters(X);
B=E_X;

%note that the second value return by FitGaussianParameters is standard
%deviation not the variance, so in the situation where we need to involve
%the variance for computation, we should transfer the deviation to variance
for i=2:N+1
   [miu,V_XUi]=FitGaussianParameters(X+U(:,i-1));
   E_XUi=(V_XUi^2-V_U(i-1)^2-V_X^2)/2+E_U(i-1)*E_X;
   B=[B;E_XUi];
end

% solve A*Beta = B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta=A\B;
% then compute sigma according to eq. (11) in PA description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE

sum_betaij_covUij=cov_Uij;
for i=1:N
    sum_betaij_covUij(i,:)=sum_betaij_covUij(i,:)*Beta(i);
end
for j=1:N
   sum_betaij_covUij(:,j)=sum_betaij_covUij(:,j)*Beta(j);
end

sum_betaij_covUij=sum(sum(sum_betaij_covUij));

%V_X is a standard deviation not variance 
sigma=sqrt(V_X^2-sum_betaij_covUij);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%