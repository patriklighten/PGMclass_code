% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeMEU( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  
  % We assume I has a single decision node.
  % You may assume that there is a unique optimal decision.
  D = I.DecisionFactors(1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  % 
  % Some other information that might be useful for some implementations
  % (note that there are multiple ways to implement this):
  % 1.  It is probably easiest to think of two cases - D has parents and D 
  %     has no parents.
  % 2.  You may find the Matlab/Octave function setdiff useful.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    EUF = CalculateExpectedUtilityFactor( I );
    OptimalDecisionRule=I.DecisionFactors;
    OptimalDecisionRule.val=OptimalDecisionRule.val*0;
    if(length(EUF.var)==1)
        [MEU,index]=max(EUF.val);
        OptimalDecisionRule.val=OptimalDecisionRule.val*0;
        OptimalDecisionRule.val(index)=1;
    else 
        [B.var, mapB] = setdiff(EUF.var, D.var(1));
        [dummy, mapD]=ismember(D.var(1),EUF.var);
         B.card =EUF.card(mapB);
         B.val = zeros(1, prod(B.card));
         assignments = IndexToAssignment(1:length(EUF.val), EUF.card);
         indxB = AssignmentToIndex(assignments(:, mapB), B.card);
         N=prod(B.card);
         
         optrules=zeros(1,N);
         MEU=ones(1,N)*(-1)+min(EUF.val);
         for i=1:length(assignments)
             if(EUF.val(i)>MEU(indxB(i)))
                 MEU(indxB(i))=EUF.val(i);
                 
                OptimalDecisionRule= SetValueOfAssignment(OptimalDecisionRule,assignments(i,:) ,1);
               
                %record last the max decistion assigment 
                %set the last max assignment rule as 0
                if(optrules(indxB(i))~=0)
                assignment1=assignments(i,:);
                assignment1(mapD)=optrules(indxB(i));
                 OptimalDecisionRule= SetValueOfAssignment(OptimalDecisionRule,assignment1 ,0);
                 
                end
                optrules(indxB(i))=assignments(i,mapD);
                 
             end 
         end
         MEU=sum(MEU);
    end
    
         
        
        
            
        
        
    end
    
    
   
    
    


