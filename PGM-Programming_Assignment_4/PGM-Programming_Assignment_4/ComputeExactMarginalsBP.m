%COMPUTEEXACTMARGINALSBP Runs exact inference and returns the marginals
%over all the variables (if isMax == 0) or the max-marginals (if isMax == 1). 
%
%   M = COMPUTEEXACTMARGINALSBP(F, E, isMax) takes a list of factors F,
%   evidence E, and a flag isMax, runs exact inference and returns the
%   final marginals for the variables in the network. If isMax is 1, then
%   it runs exact MAP inference, otherwise exact inference (sum-prod).
%   It returns an array of size equal to the number of variables in the 
%   network where M(i) represents the ith variable and M(i).val represents 
%   the marginals of the ith variable. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function M = ComputeExactMarginalsBP(F, E, isMax)

% initialization
% you should set it to the correct value in your code
M = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Implement Exact and MAP Inference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isMax==0)
    
   count=length(unique([F.var])); 
   M=repmat(struct('var', [], 'card', [], 'val', []), count, 1);
   
   count=0;
   P = CreateCliqueTree(F, E);
   [P,messages] = CliqueTreeCalibrate(P, 0);
     N=length(P.cliqueList);
     for i=1:N
         if(i==5)
             i=5;
         end
         
         for j=1:length(P.cliqueList(i).var)
             M(P.cliqueList(i).var(j))=FactorMarginalization(P.cliqueList(i),setdiff(P.cliqueList(i).var,P.cliqueList(i).var(j)));
             M(P.cliqueList(i).var(j)).val=M(P.cliqueList(i).var(j)).val/sum(M(P.cliqueList(i).var(j)));
         end
     end
end

if(isMax==1)
     
      P = CreateCliqueTree(F, E);
      
      count=length(unique([F.var])); 
     M=repmat(struct('var', [], 'card', [], 'val', []), count, 1);
      
     [P,messages] = CliqueTreeCalibrate(P, 1);
     N=length(P.cliqueList);
     for i=1:N
         for j=1:length(P.cliqueList(i).var)
             M(P.cliqueList(i).var(j))=FactorMaxMarginalization(P.cliqueList(i),setdiff(P.cliqueList(i).var,P.cliqueList(i).var(j)));
         end
     end
     
     
    
    
    
end



end
