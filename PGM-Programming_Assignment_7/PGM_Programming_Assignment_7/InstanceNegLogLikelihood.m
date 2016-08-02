% function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)
% returns the negative log-likelihood and its gradient, given a CRF with parameters theta,
% on data (X, y).                       
%
% Inputs:
% X            Data.                           (numCharacters x numImageFeatures matrix)
%              X(:,1) is all ones, i.e., it encodes the intercept/bias term.
% y            Data labels.                    (numCharacters x 1 vector)
% theta        CRF weights/parameters.         (numParams x 1 vector)
%              These are shared among the various singleton / pairwise features.
% modelParams  Struct with three fields:
%   .numHiddenStates     in our case, set to 26 (26 possible characters)
%   .numObservedStates   in our case, set to 2  (each pixel is either on or off)
%   .lambda              the regularization parameter lambda
%
% Outputs:
% nll          Negative log-likelihood of the data.    (scalar)
% grad         Gradient of nll with respect to theta   (numParams x 1 vector)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)

    % featureSet is a struct with two fields:
    %    .numParams - the number of parameters in the CRF (this is not numImageFeatures
    %                 nor numFeatures, because of parameter sharing)
    %    .features  - an array comprising the features in the CRF.
    %
    % Each feature is a binary indicator variable, represented by a struct 
    % with three fields:
    %    .var          - a vector containing the variables in the scope of this feature
    %    .assignment   - the assignment that this indicator variable corresponds to
    %    .paramIdx     - the index in theta that this feature corresponds to
    %
    % For example, if we have:
    %   
    %   feature = struct('var', [2 3], 'assignment', [5 6], 'paramIdx', 8);
    %
    % then feature is an indicator function over X_2 and X_3, which takes on a value of 1
    % if X_2 = 5 and X_3 = 6 (which would be 'e' and 'f'), and 0 otherwise. 
    % Its contribution to the log-likelihood would be theta(8) if it's 1, and 0 otherwise.
    %
    % If you're interested in the implementation details of CRFs, 
    % feel free to read through GenerateAllFeatures.m and the functions it calls!
    % For the purposes of this assignment, though, you don't
    % have to understand how this code works. (It's complicated.)
    
    featureSet = GenerateAllFeatures(X, modelParams);

    % Use the featureSet to calculate nll and grad.
    % This is the main part of the assignment, and it is very tricky - be careful!
    % You might want to code up your own numerical gradient checker to make sure
    % your answers are correct.
    %
    % Hint: you can use CliqueTreeCalibrate to calculate logZ effectively. 
    %       We have halfway-modified CliqueTreeCalibrate; complete our implementation 
    %       if you want to use it to compute logZ.
    
    nll = 0;
    grad = zeros(size(theta));
    %%%
    % Your code here:
    numhidden_state=modelParams.numHiddenStates;
    
    % there are three kinds of feature, but we only count the unique var
    % field for converting the featureset to factor list 
    numfact=size(X,1)+size(X,1)-1; %correspond for one variable features and two vairable features (Yi,Yi+1)
    factor_list=repmat(struct('var', [], 'card', [], 'val', []), 1,numfact);
    
    %we assume that the first size(X,1) facors in factor list has single
    %variable, and last size(X,1)-1 has two variables, so, for a feature is
    %it var field is [a], then it should be load in the ath factor, if
    %its var field is [a,a+1], then it should be in the size(X,1)+a th
    %factor 
    for i=1:numfact
        if(i<size(X,1)+1)
            factor_list(i).var=i;
            factor_list(i).card=modelParams.numHiddenStates;
            factor_list(i).val=zeros(1,modelParams.numHiddenStates);
        else 
            factor_list(i).var=[i-size(X,1),i+1-size(X,1)];
            factor_list(i).card=[1,1]*modelParams.numHiddenStates;
            factor_list(i).val=zeros(1,prod(factor_list(i).card));
        end
    end
    
    
    for i=1:length(featureSet.features)
        feature=featureSet.features(i);
        if (length(feature.var)==1)
            if(factor_list(feature.var).val(feature.assignment)==0)
                factor_list(feature.var).val(feature.assignment)=exp(theta(feature.paramIdx));
            else 
                factor_list(feature.var).val(feature.assignment)=factor_list(feature.var).val(feature.assignment)*exp(theta(feature.paramIdx));
            end
        else  % the var field is[x,x+1]     
            index=size(X,1)+feature.var(1);
            index1=AssignmentToIndex(feature.assignment, factor_list(index).card);
            if(factor_list(index).val(index1)==0)
               factor_list(index).val(index1)=exp(theta(feature.paramIdx));
            else
                factor_list(index).val(index1)=factor_list(index).val(index1)*exp(theta(feature.paramIdx));
            end
        end
        
                
            
    end
    
    
    %after connvert the featureset as factor list, we can use the factor to
    %build clique tree and calibrate it, to also compute the logZ 
    P = CreateCliqueTree(factor_list);
    [P, logZ] = CliqueTreeCalibrate(P, 0);
    
    %compute fature count, weight feature count, model feature count 
       featcount=zeros(1,featureSet.numParams);
       
   % there are three kinds of features, so three kinds of the parameters
   % need to be taken into account 
   % the first kind is conditional singleton feature, which has the number
   % of numCharacters * numImageFeatures * numHiddenStates;  the second
   % kind is unconditional singleton feature which has the number of
   % numCharacters*numHiddenStates;the third kind is unconditional pair
   % wise feature which has number of numCharacters * numImageFeatures^2
   
   %according to the GenerateAllFeatures, the first part of the generated
   %feature part is conditional singleton feature, then the unconditional
   %singleton feature, and finally the unconditional pair-wise feature
       
    for i=1:length(featureSet.features)
        
        %size(X,1)*numhidden_state+numhidden_state*numhidden_state*(size(X,1)-1)
        
            if(sum(featureSet.features(i).assignment==y(featureSet.features(i).var))==length(featureSet.features(i).var))
            featcount(featureSet.features(i).paramIdx)=featcount(featureSet.features(i).paramIdx)+1;
            end
            
    end
    
    nll= logZ-sum(featcount.*theta )+ modelParams.lambda/2*sum(theta.^2);
    
    
    %compute all the marginal
    M=computeallmarginal(P,factor_list,0);
    
    for i=1:length(featureSet.features)
        
       if(length(featureSet.features(i).var)==2)
           PY_x=P.cliqueList(featureSet.features(i).var(1)).val(AssignmentToIndex(featureSet.features(i).assignment, [numhidden_state,numhidden_state]));
     
       else 
           %if length(featureSet.features(i).var)==1
          PY_x=M(featureSet.features(i).var).val(featureSet.features(i).assignment);
       end
       index=featureSet.features(i).paramIdx;
       grad(index)= grad(featureSet.features(i).paramIdx)+PY_x;
    end
    
    for i=1:length(grad)
        grad(i)= grad(i)-featcount(i)+modelParams.lambda*theta(i);
    end
    
            
    
    
    
    
    
    
    
    
    

       end

       
       function M=computeallmarginal(cliquetree,factor_list,isMax)
       F=factor_list;
       PCalibrated = cliquetree;
varsList = unique([F.var]);
M = repmat(struct('var', 0, 'card', 0, 'val', []), length(varsList), 1);
for i = 1:length(varsList)
    % Iterate through variables and find the marginal for each
    clique = struct('var', 0, 'card', 0, 'val', []);
    currentVar = varsList(i);
    for j = 1:length(PCalibrated.cliqueList)
        % Find a clique with the variable of interest
        if ~isempty(find(ismember(PCalibrated.cliqueList(j).var, currentVar)))
            % A clique with the variable has been indentified
            clique = PCalibrated.cliqueList(j);
            break
        end
    end
    if isMax == 0
        % Do sum-product inference
        M(i) = FactorMarginalization(clique, setdiff(clique.var, currentVar));
        if any(M(i).val ~= 0)
            % Normalize
            M(i).val = M(i).val/sum(M(i).val);
        end
    else
        % Do MAP inference by using FactorMaxMarginalization instead of
        % FactorMarginalization
        M(i) = FactorMaxMarginalization(clique, setdiff(clique.var, currentVar));
    end
end
       end 
