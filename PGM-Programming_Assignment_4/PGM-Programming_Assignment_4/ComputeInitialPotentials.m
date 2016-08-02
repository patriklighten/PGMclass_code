%COMPUTEINITIALPOTENTIALS Sets up the cliques in the clique tree that is
%passed in as a parameter.
%
%   P = COMPUTEINITIALPOTENTIALS(C) Takes the clique tree skeleton C which is a
%   struct with three fields:
%   - nodes: cell array representing the cliques in the tree.
%   - edges: represents the adjacency matrix of the tree.
%   - factorList: represents the list of factors that were used to build
%   the tree. 
%   
%   It returns the standard form of a clique tree P that we will use through 
%   the rest of the assigment. P is struct with two fields:
%   - cliqueList: represents an array of cliques with appropriate factors 
%   from factorList assigned to each clique. Where the .val of each clique
%   is initialized to the initial potential of that clique.
%   - edges: represents the adjacency matrix of the tree. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function P = ComputeInitialPotentials(C)

% number of cliques
N = length(C.nodes);

% initialize cluster potentials 
P.cliqueList = repmat(struct('var', [], 'card', [], 'val', []), N, 1);
P.edges = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% First, compute an assignment of factors from factorList to cliques. 
% Then use that assignment to initialize the cliques in cliqueList to 
% their initial potentials. 

% C.nodes is a list of cliques.
% So in your code, you should start with: P.cliqueList(i).var = C.nodes{i};
% Print out C to get a better understanding of its structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%record the index of clique each factors are assigned to, note that  a
%factor can assign to only one clique and be assigned only one time. 
factor_toClique=zeros(1,length(C.factorList));

for i=1:N
P.cliqueList(i).var=C.nodes{i};
newFactor = struct('var', [], 'card', [], 'val', []);
if (i==9)
    newFactor=newFactor;
end
%get the joint distribution of the variables in the scope of the cluster 
for j=1:length(C.factorList)
   
    if(sum(ismember(C.factorList(j).var,C.nodes{i}))== length(C.factorList(j).var))
        %check if the factor has been assigned to a clique, if it is, don't
        %assign the factor to the clique any longer 
        if(factor_toClique(j)==0)
            factor_toClique(j)=i;
            newFactor=FactorProduct(newFactor,C.factorList(j));
        end

    end
    
end
% because the newFactor may re-arrange the variables, so assign the
% newFactor to the P.cliqueList(i)
P.cliqueList(i)=newFactor;
end

%because clusters index in P is same with C, so we can directly copy the
%edges of C to P

P.edges=C.edges;

end

