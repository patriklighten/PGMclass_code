function [P loglikelihood] = LearnCPDsGivenGraph(dataset, G, labels)
%
% Inputs:
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% G: graph parameterization as explained in PA description
% labels: N x 2 true class labels for the examples. labels(i,j)=1 if the 
%         the ith example belongs to class j and 0 elsewhere        
%
% Outputs:
% P: struct array parameters (explained in PA description)
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
K = size(labels,2);

loglikelihood = 0;
P.c = zeros(1,K);

% estimate parameters
% fill in P.c, MLE for class probabilities
% fill in P.clg for each body part and each class
% choose the right parameterization based on G(i,1)
% compute the likelihood - you may want to use ComputeLogLikelihood.m
% you just implemented.
%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
[mu sigma] = FitGaussianParameters(labels);
P.c=mu;

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
P.clg=repmat(struct('mu_x',[],'mu_y',[],'mu_angle',[],'theta',[],'sigma_x',[],'sigma_y',[],'sigma_angle',[]),size(G,1));
k=K;


for i=1:size(G,1)
    flag=0;
    
            
    
    for j=1:K
        
        if(length(size(G))==2)
          if(G(i,1)~=0)
            parentindex=G(i,2);
           else
            flag=1;
           end
        else 
            if(G(i,j,1)~=0)
               parentindex=G(i,j,2);
            else
                flag=1;
            end
        end
    dataindx=find(labels(:,j)==1);
    
      if(flag==0)
        datak=squeeze(dataset(dataindx,parentindex,:));
        theta=[];
        sigma_y=[];
        
        [t1,t2]=FitLinearGaussianParameters(dataset(dataindx,i,1),datak );
        t1=t1';
        t1=t1([length(t1),1:length(t1)-1]);
        theta=[theta,t1];
        P.clg(i).sigma_y=[P.clg(i).sigma_y,t2];
        
        [t1,t2]=FitLinearGaussianParameters(dataset(dataindx,i,2),datak );
         t1=t1';
         %note that the arrangement of t1 is different with t1 computed
         %above, the t1(N+1) is the t1(0) in in the clg, see the
         %homework description for detail
        t1=t1([length(t1),1:length(t1)-1]);
        theta=[theta,t1];
        P.clg(i).sigma_x=[P.clg(i).sigma_x,t2];
        
        [t1,t2]=FitLinearGaussianParameters(dataset(dataindx,i,3),datak );
        t1=t1';
        t1=t1([length(t1),1:length(t1)-1]);
        theta=[theta,t1];
        P.clg(i).sigma_angle=[P.clg(i).sigma_angle,t2];
        if(any(size(P.clg(i).theta))==0)
            P.clg(i).theta=theta;
        else 
        P.clg(i).theta=[P.clg(i).theta;theta];
        end
        
        
    else
        [t1,t2]=FitGaussianParameters(dataset(dataindx,i,1));
        P.clg(i).mu_y=[P.clg(i).mu_y,t1];
        P.clg(i).sigma_y=[P.clg(i).sigma_y,t2];
        [t1,t2]=FitGaussianParameters(dataset(dataindx,i,2));
        P.clg(i).mu_x=[P.clg(i).mu_x,t1];
        P.clg(i).sigma_x=[P.clg(i).sigma_x,t2];
        [t1,t2]=FitGaussianParameters(dataset(dataindx,i,3));
        P.clg(i).mu_angle=[P.clg(i).mu_angle,t1];
        P.clg(i).sigma_angle=[P.clg(i).sigma_angle,t2];
    end
    
   end
    
 end
    


 loglikelihood = ComputeLogLikelihood(P, G, dataset);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('log likelihood: %f\n', loglikelihood);

end
