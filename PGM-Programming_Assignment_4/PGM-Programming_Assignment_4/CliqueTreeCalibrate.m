%CLIQUETREECALIBRATE Performs sum-product or max-product algorithm for 
%clique tree calibration.

%   P = CLIQUETREECALIBRATE(P, isMax) calibrates a given clique tree, P 
%   according to the value of isMax flag. If isMax is 1, it uses max-sum
%   message passing, otherwise uses sum-product. This function 
%   returns the clique tree where the .val for each clique in .cliqueList
%   is set to the final calibrated potentials.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function [P,messages] = CliqueTreeCalibrate(P, isMax)


% Number of cliques in the tree.
N = length(P.cliqueList);

% Setting up the messages that will be passed.
% MESSAGES(i,j) represents the message going from clique i to clique j. 
MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We have split the coding part for this function in two chunks with
% specific comments. This will make implementation much easier.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
% While there are ready cliques to pass messages between, keep passing
% messages. Use GetNextCliques to find cliques to pass messages between.
% Once you have clique i that is ready to send message to clique
% j, compute the message and put it in MESSAGES(i,j).
% Remember that you only need an upward pass and a downward pass.
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %downstream message passing and initializing the message, all the message
 %is initialize as 1 
 
 
 
 messages=MESSAGES;
      
  

 for upstream=1:2*(N-1)
    
[i, j] = GetNextCliques(P, messages);

if(i==0)
        break;
    end
%i=downward_order(N-upstream,2);
%j=downward_order(N-upstream,1);


%initial beta_i is equal to phi_i;
beta_i=P.cliqueList(i);
edges_i= find(P.edges(i,:)==1);
if(isMax==0)
 for j1=1:sum(P.edges(i,:))
    %normalize the messages before it passed 
    if(edges_i(j1)~=j&&size(messages(edges_i(j1),i).var,2)~=0)
    messages(edges_i(j1),i).val=messages(edges_i(j1),i).val/sum(messages(edges_i(j1),i).val);
    beta_i=FactorProduct(beta_i, messages(edges_i(j1),i));
    
    end
 end
message= FactorMarginalization(beta_i,setdiff(P.cliqueList(i).var, intersect(P.cliqueList(i).var,P.cliqueList(j).var)));
message.val=message.val/sum(message.val);
messages(i,j)=message;

elseif(isMax==1) 
    %take logrithm   
    beta_i.val=log(beta_i.val);
    for j1=1:sum(P.edges(i,:))
        
    if(edges_i(j1)~=j&&size(messages(edges_i(j1),i).var,2)~=0)
    %replace the factor product as log sum factor product 
    beta_i=LogFactorProduct(beta_i,messages(edges_i(j1),i));
    
    end
    end
 
    message= FactorMaxMarginalization(beta_i,setdiff(P.cliqueList(i).var, intersect(P.cliqueList(i).var,P.cliqueList(j).var)));
    messages(i,j)=message;
end 
 end
 

 
 
MESSAGES=messages;


 
 
     
     
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Now the clique tree has been calibrated. 
% Compute the final potentials for the cliques and place them in P.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isMax==0)
for i1=1:length(P.cliqueList)
    edges=find(P.edges(i1,:)==1);
    
    for j1=1:sum(P.edges(i1,:))
        if(size(messages(edges(j1),i1).var,2)~=0)
        messages(edges(j1),i1).val=messages(edges(j1),i1).val/sum(messages(edges(j1),i1).val);
        P.cliqueList(i1)=FactorProduct(P.cliqueList(i1), messages(edges(j1),i1));
        end
    end
end

elseif(isMax==1)
 for i1=1:length(P.cliqueList)
    edges=find(P.edges(i1,:)==1);
    %take logrithm of the val field  
    P.cliqueList(i1).val=log(P.cliqueList(i1).val);
    for j1=1:sum(P.edges(i1,:))
        if(size(messages(edges(j1),i1).var,2)~=0)
            %replace the factor product as LogFactorProduct, so the 
            % final believe as well as the message are in log-space   
        P.cliqueList(i1)=LogFactorProduct(P.cliqueList(i1),messages(edges(j1),i1));
        end
    end
 end
    
end 
return

end

function [downstream_mark,downstream_node, MESSAGES,newcliques,count,node_degree,downward_order]=compute_message(P,downstream_mark,downstream_node, MESSAGES, count,node_degree,downward_order, index)
        
      if (index==9)
          index=index;
      end 
        edges=find(P.edges(index,:)==1); 
        beta_i=P.cliqueList(index);
        newcliques=[];
         % the clique that the node(index) will send message to. 
        sendcliques=[];
        for i=1:length(edges)
            if(downstream_mark(index,edges(i))==0)
               sendcliques=[sendcliques,edges(i)];
            end 
            
        end 
        
        for j1=1:length(sendcliques)
            j=sendcliques(j1);
            for i1=1:length(edges)
                if( edges(i1)~=j)
                    if(size(MESSAGES(edges(i1),index).var,2)~=0)
                        MESSAGES(edges(i1),index).val=MESSAGES(edges(i1),index).val/sum(MESSAGES(edges(i1),index).val);
                        beta_i=FactorProduct(beta_i,MESSAGES(edges(i1),index));
                  %  else 
                  %      message=FactorMarginalization(P.cliqueList(edges(i1)),setdiff(P.cliqueList(edges(i1)).var, intersect(P.cliqueList(index).var,P.cliqueList(edges(i1)).var)));
                  %      message.val=message.val/sum(message.val);                 
                  %      beta_i=FactorProduct(beta_i,message);
                    end 
                
                    
                end
            end
            MESSAGES(index,j)=FactorMarginalization(beta_i,setdiff(P.cliqueList(index).var, intersect(P.cliqueList(index).var,P.cliqueList(j).var)));
     
            %MESSAGES(index,j).val=ones(1,prod(MESSAGES(index,j).card));
            count=count+1;
            %record downward messages ordering 
            downward_order(count,:)=[index,j];
            
            downstream_mark(index,j)=1;
            downstream_mark(j,index)=1;
            
            node_degree(index)=node_degree(index)-1;
            node_degree(j)=node_degree(j)-1;
        end 
        
        for i=1:length(edges)
            if(node_degree(edges(i))>0)
                newcliques=[newcliques,edges(i)];
            else 
                downstream_node(edges(i))=1;
            end
            
        end 
      
        downstream_node(index)=1;
               
        
        

end 


