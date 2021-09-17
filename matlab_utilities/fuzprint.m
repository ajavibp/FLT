% FUZPRINT Converts a fuzzy model in its linguistic representation. 
%
% Use the linguistic form:
%   "IF Name_of_Input1 is...and Name_of_Input2 is...THEN Name_of_Output1 is..."
% 
%   fuzprint('File.txt',Fuzzy_Model)
% 
%   fuzprint('File.txt',Fuzzy_Model,Name_of_Input1,Name_of_Input2,...,
%            Name_of_Output1,Name_of_Output2,...)
% 
%   fuzprint('File.txt',Fuzzy_Model,Name_of_Input1,Name_of_Input2,...,
%            Name_of_Output1,Name_of_Output2,...,Accuracy)
% 
% Arguments:
% 
%   'File.txt' -> String with the name of file to write. This file will be
%                   overwrited without confirmation.
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
%   Name_of_InputX -> String with the name of input X.
%
%   Name_of_OutputX -> String with the name of output X.
%
%   Accuracy -> Number of decimals to print (10 by defult).
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzeval, fuzjac, fuzlinear, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
