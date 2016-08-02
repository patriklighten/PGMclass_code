%GETNEXTCLIQUES Find a pair of cliques ready for message passing
%   [i, j] = GETNEXTCLIQUES(P, messages) finds ready cliques in a given
%   clique tree, P, and a matrix of current messages. Returns indices i and j
%   such that clique i is ready to transmit a message to clique j.
%
%   We are doing clique tree message passing, so
%   do not return (i,j) if clique i has already passed a message to clique j.
%
%	 messages is a n x n matrix of passed messages, where messages(i,j)
% 	 represents the message going from clique i to clique j. 
%   This matrix is initialized in CliqueTreeCalibrate as such:
%      MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);
%
%   If more than one message is ready to be transmitted, return 
%   the pair (i,j) that is numerically smallest. If you use an outer
%   for loop over i and an inner for loop over j, breaking when you find a 
%   ready pair of cliques, you will get the right answer.
%
%   If no such cliques exist, returns i = j = 0.
%
%   See also CLIQUETREECALIBRATE
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function [i, j] = GetNextCliques(P, messages)

% initialization
% you should set them to the correct values in your code
i = 0;
j = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(messages,1);

message_set=zeros(N,2);

count=0;
count1=0;
J=0;

for i1=1:N
    count=0;
    %this function only return one message transmitting path 
    if(count1>0)
        break;
    end
    if(sum(P.edges(i1,:))>0)
        edges= find(P.edges(i1,:)==1);
        for j1=1:sum(P.edges(i1,:))
            %the if statement check the message from the clique's downstream neighbors 
            if(size(messages(edges(j1),i1).var,2)==0)
                count=count+1;
                J=edges(j1);
            end
        end
       if(count==1&&size(messages(i1,J).var,2)==0)
           count1=count1+1;
           message_set(count1,:)=[i1,J];
       elseif(count==0)
           %if count==0, it means the clique has recieved all the neighbors
           %in downstream and now the clique need to send massage to one of
           %its neighbor first which the clique has no message to it. 
           for j1=1:sum(P.edges(i1,:))
               
                if(size(messages(i1,edges(j1)).var,2)==0)
                    count1=count1+1;
                    message_set(count1,:)=[i1,edges(j1)];
                    
                    break;
                end
           end
       end
       
    end
    
end

 %If more than one message is ready to be transmitted, return 
%   the pair (i,j) that is numerically smallest
if (count1>0)
    i=message_set(1,1);
    j=message_set(1,2);
end


  
    
           
       


return;
