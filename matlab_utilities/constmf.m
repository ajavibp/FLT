function y=constmf(x,c)
% CONSTMF Constant membership function.
%
%   constmf(x,c) always returns c.
%
% For example:
%
%   x = (0:0.1:10)';
%   y = constmf(x, 0.5);
%   plot(x, y);
%
% See also dsigmf, evalmf, gauss2mf, gaussmf, gbellmf, mf2mf, anymf, pimf
%          psigmf, sigmf, smf, trimf, trapmf, zmf.

if nargin ~= 2
    error('Two arguments are required by the Const membership function.');
elseif length(C) ~= 1
    error('2nd parameter of the function must be a scalar.');
end
y = ones(size(x))*C;