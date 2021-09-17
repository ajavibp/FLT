% FUZLINEAR Calculartes the linealization of a Fuzzy Model in a point.
%
%   A = fuzlinear(Fuzzy_Model, X, U)
%
%   [A,B] = fuzlinear(Fuzzy_Model, X, U)
%
%   [A,B,F] = fuzlinear(Fuzzy_Model, X, U)
%
% Solve the matrices the linear representation of the 'Fuzzy_Model' model
% according to Taylor series in the point X with signal control U.
% F is the value of the function at the point (X,U).
%
%                           Y = F + AX + BU
%
% Arguments:
%
%   Fuzzy_Model -> Input Fuzzy model.
%
%   X -> State vector.
%
%   U -> Control vector.
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat, fuzcomb, fuzeval, fuzjac, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
