% FUZJAC Calculates the Jacobian matrix of the closed loop fuzzy system.
%
%   J = fuzjac(X, Plant, Controller)
%
% Arguments:
%
%   X -> State vector where the jacobian matrix is calculated.
%
%   Plant -> Fuzzy model of the Plant.
%
%   Controlador -> Fuzzy model of the Controller.
%
%   J -> Jacobian matrix of closed loop system.
%
%                                  (  df1          df1  )
%     dx1                          | -----   ...  ----- |
%    ----- = f1(x)                 |  dx1          dxn  |
%     dt                           |                    |
%     ...                      J = |  ...    ...   ...  |
%     dxn                          |                    |
%    ----- = fn(x)                 |  dfn          dfn  |
%     dt                           | -----   ...  ----- |
%                                  (  dx1          dxn  )
%    x = [ x1 ; x2 ; ... ; xn]
%
%   Fuzzy_Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%                  or a 'FIS' variable from MATLAB Workspace.
%
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%         fuz2mat, fuzcomb, fuzeval, fuzlinear, fuzprint, mat2antec, mat2conseq,
%         mat2fuz, fuzsubsystem, txt2fis
