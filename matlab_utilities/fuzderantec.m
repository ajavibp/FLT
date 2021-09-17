% FUZDERANTEC Gets the derivatives of a fuzzy model respect its antecedents.
%
%   It is assumed that the output of the system is the same as the parameter,
%   otherwise it should be considered derivative equal to 0.
%
%   dh_dantec = fuzderconseq(Fuzzy_Model, Point)
%
%   dh_dantec = fuzderconseq(Fuzzy_Model, Point, input, out, rule, parameter)
%
% Arguments:
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
%   Point -> Point where the derivative will be calculate.
%
%   input -> Input for the parameter.
%
%   out -> Output for select the rule for the parameter.
%
%   rule -> Rule for the output 'out' for the parameter.
%
%   parameter -> Parameter of the antecedent.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzderconseq, fuzderparam, fuzeval, fuzjac,
%          fuzlinear, fuzprint, mat2antec, mat2conseq, mat2fuz, fuzsubsystem,
%          txt2fis
