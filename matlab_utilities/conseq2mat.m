% CONSEQ2MAT Extracts data from consequent of a fuzzy model.
% 
%   Vector = conseq2mat(Fuzzy_Model)
% 
%   Vector = conseq2mat(Fuzzy_Model,Output,Rule)
% 
%   [Vector, length] = conseq2mat(Fuzzy_Model)
%
%   [Vector, length] = conseq2mat(Fuzzy_Model,Output,Rule)
%
% Arguments:
% 
%   Vector -> Vector with data from fuzzy model consequent.
%
%   Output -> Use to select the consequent of one output.
%
%   Rule -> Use to select the consequent of one rule for given output.
%
%   length -> The length of Vector, 0 if error.
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%     
% See also activation, antec2mat, aproxjac, aproxlinear, fis2txt, fuz2mat,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
