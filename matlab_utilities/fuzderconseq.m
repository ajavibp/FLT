% FUZDERCONSEQ Gets the derivative of a fuzzy model respect its consequents.
%
%   It is assumed that the output of the system is the same as the parameter,
%   otherwise it should be considered derivative equal to 0.
%
%   dh_dconseq = fuzderconseq(Fuzzy_Model, Point)
%
%   dh_dconseq = fuzderconseq(Fuzzy_Model, Point, out, rule, parameter)
%
% Arguments:
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
%   Point -> Point where the derivative will be calculate.
%
%   out -> Output for select the rule for the parameter.
%
%   rule -> Rule for the output 'out' for the parameter.
%
%   parameter -> Parameter of the consequent (0 -> affine consequent)
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzderantec, fuzderparam, fuzeval, fuzjac, fuzlinear,
%          fuzprint, mat2antec, mat2conseq, mat2fuz, fuzsubsystem, txt2fis
