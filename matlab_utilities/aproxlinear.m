% APROXLINEAR Linealizes a Fuzzy Model usign finite differences.
%
%   A = aproxlinear(Fuzzy_Model, X, U)
%
%   A = aproxlinear(Fuzzy_Model, X, U, h)
%
%   [A,B] = aproxlinear(Fuzzy_Model, X, U, h)
%
%   [A,B,F] = aproxlinear(Fuzzy_Model, X, U, h)
%
% Solve the matrices the linear representation of the 'Fuzzy_Model' model
% usign finite differences in the point X with signal control U.
% F is the value of the function at the point (X,U).
%
%                           Y = F + Ax + Bu
%
% Arguments:
%
%   Fuzzy_Model -> Input Fuzzy model.
%
%   X -> State vector.
%
%   U -> Control vector.
%
%   h -> Increment. Default value is 0.001;
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
% See also activation, antec2mat, aproxjac, conseq2mat, fis2txt, fuz2mat,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
