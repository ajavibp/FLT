% FUZCOMB Combines some fuzzy models in one.
% 
%   fuzcomb(Combined_Model_TXT,Model1, Model2,...)
%
%   Combined_Model_FIS = fuzcomb(Model1, Model2,...)
%
% Arguments:
% 
%   Combined_Model_TXT -> String with the name of TXT file to write with the
%                         combined model. This file will be overwrited without
%                         confirmation.
%
%   Combined_Model_FIS -> Variable 'FIS' to write the combined model.
%
%   ModelX -> Fuzzy models to combine. These models must have the same number
%             of inputs. They could be a '.txt' or '.fis' file, or a 'FIS'
%             variable from MATLAB Workspace.
%
% To simulate this combined model should use 'fuzeval' function, 'evalfis'
% doesn't work correctly.
%     
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
