% FUZfuzsubsystem Generates a fuzsubsystem (submodel) from a fuzzy model.
% 
%   SubFIS = fuzfuzsubsystem(FIS_Model,Outputs,Rules)
% 
% Arguments:
%
%   FIS_Model -> Original fuzzy model.
%
%   Outputs -> List of outputs for selected rules.
%
%   Rules -> Selected rules (for each output) to make the fuzsubsystem.
%
%   SubFIS -> Generated fuzsubsystem.
%
% For example:
%
%   P = readfis('Plant.fis')
%   subP = fuzfuzsubsystem(P,[1,2,2],[1,1,2])
%
%  'subP' has 3 rules, the 1st rule for the 1st output of P, and the 1st and
%  2nd rule for the 2nd output of P.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec,
%          mat2conseq, mat2fuz, txt2fis
