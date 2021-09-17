% MAT2CONSEQ Changes the consequent of a fuzzy model.
%
%   Final_Model = mat2conseq(Initial_Model,Vector)
%
%   Final_Model = mat2conseq(Initial_Model,Vector,Output,Rule)
%
%   [Final_Model,length] = mat2conseq(Initial_Model,Vector)
%
%   [Final_Model,length] = mat2conseq(Initial_Model,Vector,Output,Rule)
%
% Arguments:
%
%   Final_Model -> Changed fuzzy model.
%
%   Output -> Use to select the consequent of one output.
%
%   Rule -> Use to select the consequent of one rule for given output.
%
%   length -> The length of Vector, 0 if error.
%
%   Initial_Model -> Fuzzy Model to change. This model could be a '.txt' file,
%                    a '.fis' file, or a 'FIS' variable from MATLAB Workspace.
%
%   Vector -> Vector with data for fuzzy model consequent.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec,
%          mat2fuz, fuzsubsystem, txt2fis
