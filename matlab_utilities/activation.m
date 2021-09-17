% ACTIVATION Calculates the fulfillment degree and the derivative of a rule.
%
%   W = activation(Fuzzy_Model, Point, out)
%
%   [W, SumW] = activation(Fuzzy_Model, Point, out)
%
% For a specific rule:
%
%   W = activation(Fuzzy_Model, Point, out, rule)
%
%   [W, dW] = activation(Fuzzy_Model, Point, out, rule)
%
% Arguments:
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
%   out -> Output of the model for select the rule.
%
%   rule -> The rule for the output 'out'. If 'rule' is omitted will be
%             calculated the degree of activation of all rules of 'out'.
%
%   Point -> Point where the fulfillment degree will be calculate.
%
%   W -> Matching degree of the rule/s.
%
%   SumW -> Sumatory of all fulfillment degrees for the output 'out' => sum(W)
%
%   dW -> Derivative of the fulfillment degree. 'rule' must be used.
%   Note: dW will be correct only if all inputs of Fuzzy_Model are independent.
%
% See also antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt, fuz2mat,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
