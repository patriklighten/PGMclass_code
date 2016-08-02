%CREATECLUSTERGRAPH Takes in a list of factors and returns a Bethe cluster
%   graph. It also returns an assignment of factors to cliques.
%
%   C = CREATECLUSTERGRAPH(F) Takes a list of factors and creates a Bethe
%   cluster graph with nodes representing single variable clusters and
%   pairwise clusters. The value of the clusters should be initialized to 
%   the initial potential. 
%   It returns a cluster graph that has the following fields:
%   - .clusterList: a list of the cluster beliefs in this graph. These entries
%                   have the following subfields:
%     - .var:  indices of variables in the specified cluster
%     - .card: cardinality of variables in the specified cluster
%     - .val:  the cluster's beliefs about these variables
%   - .edges: A cluster adjacency matrix where edges(i,j)=1 implies clusters i
%             and j share an edge.
%  
%   NOTE: The index of the cluster for each factor should be the same within the
%   clusterList as it is within the initial list of factors.  Thus, the cluster
%   for factor F(i) should be found in P.clusterList(i) 
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function P = CreateClusterGraph(F, Evidence)
P.clusterList = [];
P.edges = [];
for j = 1:length(Evidence),
    if (Evidence(j) > 0),
        F = ObserveEvidence(F, [j, Evidence(j)]);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=0;%number of factors whose has more than one variables
variables=unique([F.var]);
for i=1:length(F)
    %firstly consider the factor that have more than one variables
    if(length(F(i).var)>1)
     if (count==0)
     P.clusterList=F(i);
     count=count+1;
     else 
         P.clusterList=[P.clusterList,F(i)];
         count=count+1;
     end
    end
    
end

%compute the card for individual vatiables 
V = unique([F(:).var]);
%V_factors=repmat(struct('var', [], 'card', [], 'val', []), 1, length(V))
% Setting up the cardinality for the variables since we only get a list 
% of factors.
V_clust=[];
for i = 1 : length(V),

	 for j = 1 : length(F)
		  if (~isempty(find(F(j).var == V(i))))
				%card(i) = F(j).card(find(F(j).var == V(i)));
                factor1=FactorMarginalization(F(j),setdiff(F(j).var,V(i)));
                factor1.val=factor1.val/sum(factor1.val);
                if(i==1)
                    V_clust=factor1;
                else
                    V_clust=[V_clust,factor1];
                end
                
				break;
		  end
	 end
end
count1=length(V);
P.clusterList=[V_clust,P.clusterList];

N=length(P.clusterList);
P.edges=zeros(N,N);
for i=1:count
    for j=1:length(variables)
         if(ismember(variables(j),P.clusterList(i+count1).var))
             P.edges(i+count1,j)=1;
             P.edges(j,i+count1)=1;
         end
    end
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

