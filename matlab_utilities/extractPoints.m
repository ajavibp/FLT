% EXTRACTPOINTS Extracts the significative points of a fuzzy model.
% 
%   Points = extractPoints(System)
%
%   Points = extractPoints(System, numMeshPoints)
% 
%   Points = extractPoints(System, numMeshPoints, precision)
% 
%   Points = extractPoints(System, numMeshPoints, precision, controlPoints)
% 
% Arguments:
% 
%   System -> This fuzzy model could be a '.txt' file, a '.fis' file,
%             or a 'FIS' variable from MATLAB Workspace.
%
%   numMeshPoints -> Used in memberships functions without significative
%                    points (ANYMF, CONSTMF, ...). Default value is 5.
%
%   precision -> Used to remove similar points. Default value is 1e-3.
%
%   controlPoints -> If controlPoints>0 extract all point, otherwise, only
%               extracts the state vector points, not the control vector points.
%               Default value is 0 (do not extract the control vector points).
%     
% See also activation, antec2mat, aproxjac, aproxlinear, conseq2mat, fis2txt,
%          fuz2mat,fuzcomb, fuzeval, fuzjac, fuzlinear, fuzprint,
%          initializeController, mat2antec, mat2conseq, mat2fuz, fuzsubsystem
%          txt2fis
