% FUZEVAL Performs fuzzy inference calculation of open/closed loop fuzzy system.
% 
%   Output = fuzeval(X_U,Fuzzy_Model)
% 
%   Output = fuzeval(X,Plant,Controller)
%
%   X_U -> Input vector: [X;U] = [x1;x2;...;xn;u1;u2;...;um]
%
%   X -> State vector.
%
%   Plant -> Fuzzy model of the Plant.
%
%   Controller -> Fuzzy model of the Controller.
%         
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
