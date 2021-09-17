% MAT2FUZ Changes all data of a fuzzy model.
%
%   Final_Model = mat2fuz(Initial_Model,Vector)
%
%   Final_Model = mat2fuz(Initial_Model,Vector,Output,Rule)
%
%   [Final_Model,length] = mat2fuz(Initial_Model,Vector)
%
%   [Final_Model,length] = mat2fuz(Initial_Model,Vector,Output,Rule)
%
% Arguments:
%
%   Final_Model -> Changed fuzzy model.
%
%   Output -> Use to select the output.
%
%   Rule -> Use to select the rule for given output.
%
%   length -> The length of Vector, 0 if error.
%
%   Initial_Model -> Fuzzy model to change. This model could be a '.txt' file,
%                    a '.fis' file, or a 'FIS' variable from MATLAB Workspace.
%
%   Vector -> Vector with data for fuzzy model.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec,
%          mat2conseq, fuzsubsystem, txt2fis
