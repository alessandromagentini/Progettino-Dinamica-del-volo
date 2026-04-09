function opt = get_integration_opt(tol, maxstep, varargin)
% Function per le opzioni di integrazione numerica

if nargin < 3
    manovra = [];
else
    manovra = varargin{1};
end
opt = odeset('RelTol', tol, 'AbsTol', tol, 'MaxStep', maxstep, 'Events', manovra);