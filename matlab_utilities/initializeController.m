% INITIALIZECONTROLLER Initializes a fuzzy Controller from a fuzzy Plant.
% 
%   Controller = initializeController(Plant)
% 
% Arguments:
% 
%   Plant -> Plant model. This model could be a '.txt' file, a '.fis' file,
%            or a 'FIS' variable from MATLAB Workspace.
%
%   Controller -> An initial controller based in the Plant. This controller
%                 has all its consecuents to 0, and each of its outputs have
%                 so many rules as all of the plant.
%     
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat,fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec,
%          mat2conseq, mat2fuz, fuzsubsystem, txt2fis
