function loglikelihood = ComputeLogLikelihood(P, G, dataset)
% returns the (natural) log-likelihood of data given the model and graph structure
%
% Inputs:
% P: struct array parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description)
%
%    NOTICE that G could be either 10x2 (same graph shared by all classes)
%    or 10x2x2 (each class has its own graph). your code should compute
%    the log-likelihood using the right graph.
%
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% 
% Output:
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset,1); % number of examples
K = length(P.c); % number of classes

loglikelihood = 0;
% You should compute the log likelihood of data as in eq. (12) and (13)
% in the PA description
% Hint: Use lognormpdf instead of log(normpdf) to prevent underflow.
%       You may use log(sum(exp(logProb))) to do addition in the original
%       space, sum(Prob).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE

%G(i; 1) = 0 indicates that body part i only has the class variable as its parent. In this
%case, G(i; 2) can take on an arbitrary value. Its parameterization follows the naive Bayes
%model described in Equations (2,3,4). The parameters can be found in the following
%12 vectors: P.clg(i).mu_x, P.clg(i).mu_y, P.clg(i).mu_angle, P.clg(i).sigma_x,
%P.clg(i).sigma_y, and P.clg(i).sigma_angle. For example, P.clg(2).sigma_y(2) is
%the standard deviation of the y-position of the head given that the class label is 2 (aliens).
% G(i; 1) = 1 indicates that body part i has, besides the class variable, another parent
%G(i; 2), and the CPDs are parameterized following Equations (5,6,7). The parameters can
%be found in the 2  12 matrix P.clg(i).theta and 1  2 vectors P.clg(i).sigma_x,
%P.clg(i).sigma_y, and P.clg(i).sigma_angle. For example, P.clg(i).theta(k, 9)
%corresponds to (9)
%ik in Equation (7)
singlelog=zeros(N,2);
for j=1:N 
     for k=1:K
         singlelog(j,k)=log(P.c(k));
        for i=1:10
            parent_data=ones(1,4);
            if(length(size(G))==2)
                G_i1=G(i,1);
                if(G_i1~=0)
                    parent_data(2:4)=dataset(j,G(i,2),:);
                end
                
            else 
                G_i1=G(i,k,1);
                if(G_i1~=0)
                    parent_data(2:4)=dataset(j,G(i,k,2),:);
                end
            end
            if(G_i1~=0)
                    mu1=sum(P.clg(i).theta(k, 1:4).*parent_data);
                    mu2=sum(P.clg(i).theta(k, 5:8).*parent_data);
                    mu3=sum(P.clg(i).theta(k, 9:12).*parent_data);
            else 
                    mu1=P.clg(i).mu_y(k);
                  
                    mu2=P.clg(i).mu_x(k);
                    mu3=P.clg(i).mu_angle(k);
            end
        singlelog(j,k)=singlelog(j,k)+lognormpdf(dataset(j,i,1),mu1,P.clg(i).sigma_y(k))+lognormpdf(dataset(j,i,2),mu2,P.clg(i).sigma_x(k))+lognormpdf(dataset(j,i,3),mu3,P.clg(i).sigma_angle(k));
                
        end 
     end
end

sumsinglelog= sum(exp(singlelog),2);
loglikelihood= sum(log(sumsinglelog));






                
                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
