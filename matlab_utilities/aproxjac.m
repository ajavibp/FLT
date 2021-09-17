% APROXJAC Calculates the closed loop Jacobian matrix usign finite differences.
%
%   J = aproxjac(X, Plant, Controller)
%
%   J = aproxjac(X, Plant, Controller, h)
%
% Arguments:
%
%   X -> State vector.
%
%   Plant -> Fuzzy model of the Plant.
%
%   Controlador -> Fuzzy model of the Controller.
%
%   h -> Increment. Default value is 0.001;
%
%   J -> Jacobian matrix of the closed loop fuzzy system.
%
%                                 (  df1          dfn  )
%    dx1                          | -----   ...  ----- |
%   ----- = f1(x)                 |  dx1          dxn  |
%    dt                           |                    |
%    ...                      J = |  ...    ...   ...  |
%    dxn                          |                    |
%   ----- = fn(x)                 |  dfn          dfn  |
%    dt                           | -----   ...  ----- |
%                                 (  dx1          dxn  )
%   x = [ x1 ; x2 ; ... ; xn]
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
% See also activation, antec2mat, aproxlinear, conseq2mat, fis2txt, fuz2mat,
%          fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint, mat2antec, mat2conseq,
%          mat2fuz, fuzsubsystem, txt2fis
