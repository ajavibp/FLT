function y = anymf(x,null)
%  ANYMF No membership function.
%     You should not use ANYMF membership funtion with MATLAB functions,
%     because MATLAB doesn't support this. Use ANYMF only with
%     'Fuzzy Logic Tools' functions.
%  
%     In 'Fuzzy Logic Tools' software, this membership function is used to
%     remove a variable from the antecedent of a rule. For example:
%
%     To make 'If x1 is A1 and x3 is A2 ...' you can use 'x2 is ANYMF'.
%
%     See also constmf, dsigmf, evalmf, gauss2mf, gaussmf, gbellmf, mf2mf, pimf
%     psigmf, sigmf, smf, trimf, trapmf, zmf.

warning('ANYMF membership function should not be executed.');

y = ones(1,length(x)); % This value is wrong, only is used to avoid errors in MATLAB functions
