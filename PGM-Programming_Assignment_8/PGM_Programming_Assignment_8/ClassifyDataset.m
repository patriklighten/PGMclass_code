function accuracy = ClassifyDataset(dataset, labels, P, G)
% returns the accuracy of the model P and graph G on the dataset 
%
% Inputs:
% dataset: N x 10 x 3, N test instances represented by 10 parts
% labels:  N x 2 true class labels for the instances.
%          labels(i,j)=1 if the ith instance belongs to class j 
% P: struct array model parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description) 
%
% Outputs:
% accuracy: fraction of correctly classified instances (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
accuracy = 0.0;

K=size(labels,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classlabel=zeros(1,N);
dataprob=zeros(1,N);
for j=1:N 
     for k=1:K
         singlelog=log(P.c(k));
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
        singlelog=singlelog+lognormpdf(dataset(j,i,1),mu1,P.clg(i).sigma_y(k))+lognormpdf(dataset(j,i,2),mu2,P.clg(i).sigma_x(k))+lognormpdf(dataset(j,i,3),mu3,P.clg(i).sigma_angle(k));
        
        end
        %after end of i =1:10 parts loop
        if(k==1)
            dataprob(j)=singlelog;
            classlabel(j)=k;
        else 
            if(dataprob(j)<singlelog)
                dataprob(j)=singlelog;
                classlabel(j)=k;
            end
        end
     end
end

count=0;
for i=1:N
    if(labels(i,classlabel(i))==1)
        count=count+1;
    end
end

accuracy=count/N;


fprintf('Accuracy: %.2f\n', accuracy);