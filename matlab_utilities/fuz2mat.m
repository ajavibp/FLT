% FUZ2MAT Extracts all data from a fuzzy model.
% 
%   Vector = fuz2mat(Fuzzy_Model)
% 
%   Vector = fuz2mat(Fuzzy_Model,Output,Rule)
% 
%   [Vector, length] = fuz2mat(Fuzzy_Model)
%
%   [Vector, length] = fuz2mat(Fuzzy_Model,Output,Rule)
%
% Arguments:
% 
%   Vector -> Vector with data from fuzzy model.
%
%   Output -> Use to select the output.
%
%   Rule -> Use to select the rule for given output.
%
%   length -> The length of Vector, 0 if error.
%
%   Fuzzy_Model -> This model could be a '.txt' file, a '.fis' file or a
%                  'FIS' variable from MATLAB Workspace.
%     
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
