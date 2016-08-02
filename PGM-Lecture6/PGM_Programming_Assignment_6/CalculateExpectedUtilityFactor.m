% Copyright (C) Daphne Koller, Stanford University, 2012

function EUF = CalculateExpectedUtilityFactor( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: A factor over the scope of the decision rule D from I that
  % gives the conditional utility given each assignment for D.var
  %
  % Note - We assume I has a single decision node and utility node.
  EUF = [];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  %
  F=[I.RandomFactors,I.UtilityFactors];
  %note that the expected utility factor's vars only include the decision node
  %and its parents nodes variables, it may not include the all variables in the utility factors  
  F1=VariableElimination(F,setdiff(unique([F.var]),[I.DecisionFactors.var]));
   
 EUF = struct('var', [], 'card', [], 'val', []);
  for i=1:length(F1)
      EUF=FactorProduct(EUF,F1(i));
  end 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


  
end  
