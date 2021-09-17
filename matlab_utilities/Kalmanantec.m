function [Model, P, err] = Kalmanantec(Model, inputs, outputs, covariance, P)
% KALMANANTEC. Calculates an iteration of EKF only for antecedents.
%
%   [Model, P, err] = Kalmanantec(Model, inputs, outputs, covariance, P)
%
%   [Model, P, err] = Kalmanantec(Model, inputs, outputs, covariance, P, Phi)
%
% Arguments:
% 
%   Model -> This fuzzy model could be a '.txt' file, a '.fis' file,
%            or a 'FIS' variable from MATLAB Workspace.
%
%   inputs -> Current input vector.
%
%   outputs -> Real outputs of the system in the current iteration.
%
%   covariance -> Noise covariance matrix estimated from the hope operator.
%
%   P -> Covariance matrix of the filter.
%
%   Phi -> Jacobian matrix that relates the parameters to be set with the
%          following value of these parameters.
%          If not specified, Phi is assumed to be the identity matrix.
%
% See also Kalmanconseq, Kalmanfuz, activation, antec2mat, aproxjac, aproxlinear
%          conseq2mat, fis2txt, fuz2mat, fuzcomb, fuzeval, fuzjac, fuzlinear,
%          fuzprint, mat2antec, mat2conseq, mat2fuz, fuzsubsystem, txt2fis
