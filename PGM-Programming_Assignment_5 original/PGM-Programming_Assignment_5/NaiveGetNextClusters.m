%NAIVEGETNEXTCLUSTERS Takes in a node adjacency matrix and returns the indices
%   of the nodes between which the m+1th message should be passed.
%
%   Output [i j]
%     i = the origin of the m+1th message
%     j = the destination of the m+1th message
%
%   This method should iterate over the messages in increasing order where
%   messages are sorted in ascending ordered by their destination index and 
%   ties are broken based on the origin index. (note: this differs from PA4's
%   ordering)
%
%   Thus, if m is 0, [i j] will be the pair of clusters with the lowest j value
%   and (of those pairs over this j) lowest i value as this is the 'first'
%   element in our ordering. (this difference is because matlab is 1-indexed)
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function [i, j] = NaiveGetNextClusters(P, m)

    i = size(P.clusterList,1);
    j = size(P.clusterList,1);
    N=sum(sum(P.edges));%total number of edges
    m=mod(m,N);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    % Find the indices between which to pass a cluster
    % The 'find' function may be useful
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %edges is an adjacency matrix
    sort_edges=zeros(sum(sum(P.edges)),2);
    count=0;
    
    
    for j1=1:size(P.clusterList,2)
        index_i=find(P.edges(:,j1)==1);
        for i1=1:sum(P.edges(:,j1))
            count=count+1;
            sort_edges(count,:)=[index_i(i1),j1];
        end
    end
    
    i=sort_edges(m+1,1);
    j=sort_edges(m+1,2);
    
end

