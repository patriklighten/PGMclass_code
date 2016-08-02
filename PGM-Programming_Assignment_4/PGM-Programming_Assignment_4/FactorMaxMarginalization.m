% FactorMaxMarginalization Max-marginalizes a factor 
% by taking the max over a given set variables.
% 
%   B = FactorMaxMarginalization(A,V) computes the factor with the variables
%   in V maxed out. The factor data structure has the following fields:
%       .var    Vector of variables in the factor, e.g. [1 2 3]
%       .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
%       .val    Value table of size prod(.card)
%
%   B.var will be A.var minus V.
%   For each assignment in B, its value is the maximum value in A 
%   of all assignments in A consistent with that assignment in B.
%
%   The resultant factor should have at least one variable remaining or this
%   function will throw an error.
%
%   This is exactly the same as FactorMarginalization, 
%   but with the sum replaced by a max.
% 
%   See also FactorMarginalization.m, IndexToAssignment.m, and AssignmentToIndex.m
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function B = FactorMaxMarginalization(A, V)

% Check for empty factor or variable list
if (isempty(A.var) || isempty(V)), B = A; return; end;

% Construct the output factor over A.var \ V (the variables in A.var that are not in V)
% and mapping between variables in A and B
[B.var, mapB] = setdiff(A.var, V);

% Check for empty resultant factor
if isempty(B.var)
  error('Error: Resultant factor has empty scope');
end;

% initialization
% you should set them to the correct values in your code
B.card = [];
B.val = [];
B.card = A.card(mapB);
assignments = IndexToAssignment(1:length(A.val), A.card);
indxB = AssignmentToIndex(assignments(:, mapB), B.card);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
% Correctly set up and populate the factor values of B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B.card = A.card(mapB);

%it should not set the initial value of B as zero, because when it is uesed
%for logrithm the log(P) always smaller than 0( which correspond to P=1)
B.val = ones(1,prod(B.card))*2;

%it should not set the initial value of maximum value for B as zero, because when it is uesed
%for logrithm the log(P) always smaller than 0 ( which correspond to P=1)
B.val = zeros(1,prod(B.card));
%so we can initialize the max_values of B as 2, so that it means it has not
%been set value yet. 
max_value=ones(1,prod(B.card))*2;
% however we can re-initialize the max_value as one of the A value as below
for i=1:length(A.val)
    if(max_value(indxB(i))==2);
        %if the max_value has not been re-initialized by A's value, we
        %initialize it 
        max_value(indxB(i))=A.val(i);
    elseif(A.val(i)>max_value(indxB(i)))
        max_value(indxB(i))=A.val(i);
    end
end


B.val=max_value;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
