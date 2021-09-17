% ANTEC2MAT Extracts data from antecedent of a fuzzy model.
% 
%   Vector = antec2mat(Fuzzy_Model)
%
%   Vector = antec2mat(Fuzzy_Model,Output,Rule)
%
%   [Vector, length] = antec2mat(Fuzzy_Model)
%
%   [Vector, length] = antec2mat(Fuzzy_Model,Output,Rule)
%
% Arguments:
% 
%   Vector -> Vector with data from fuzzy model antecedent.
%
%   Output -> Use to select the antecedent of one output.
%
%   Rule -> Use to select the antecedent of one rule for the given output.
%
%   length -> The length of Vector, 0 if error.
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%     
% See also activation, aproxjac, aproxlinear, conseq2mat, fis2txt, fuz2mat,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
