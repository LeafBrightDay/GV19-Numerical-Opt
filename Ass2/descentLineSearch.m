function [xMin, fMin, nIter, info] = descentLineSearch(F, descent, ls, alpha0, x0, tol, maxIter)
% DESCENTLINESEARCH Wrapper function executing  descent with line search
% [xMin, fMin, nIter, info] = descentLineSearch(F, descent, ls, alpha0, x0, tol, maxIter)
%
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% descent: specifies descent direction {'steepest', 'newton'}
% ls: function handle for computing the step length
% alpha0: initial step length
% rho: in (0,1) backtraking step length reduction factor
% c1: constant in sufficient decrease condition f(x_k + alpha_k*p_k) > f_k + c1*alpha_k*(df_k')*p_k)
%     Typically chosen small, (default 1e-4).
% x0: initial iterate
% tol: stopping condition on minimal allowed step
%      norm(x_k - x_k_1)/norm(x_k) < tol;
% maxIter: maximum number of iterations
%
% OUTPUTS
% xMin, fMin: minimum and value of f at the minimum
% nIter: number of iterations
% info: structure with information about the iteration
%   - xs: iterate history
%   - alphas: step lengths history
%
% Copyright (C) 2017  Marta M. Betcke, Kiko Rullan

% Parameters
% Stopping condition {'step', 'grad'}
cond1=strcmp(descent,'steepest');
cond2=strcmp(descent,'newton');
if cond1
    stopType = 'step';
end
if cond2
    stopType = 'grad';
end
% Initialization
nIter = 0;
x_k = x0;
info.xs = x0;
info.alphas = alpha0;
stopCond = false;

% Loop until convergence or maximum number of iterations
while (~stopCond && nIter <= maxIter)   
    x_k_1=info.xs(:,nIter+1);
    % =========================================================================
    switch stopType
        case 'step'
            p_k_1=-F.df(x_k_1)./norm(F.df(x_k_1));
            alpha=ls(x_k_1,p_k_1,info.alphas(nIter+1));
            x_k=x_k_1+alpha.*p_k_1;
            % Compute relative step length
            normStep = norm(x_k - x_k_1)/norm(x_k_1);
            stopCond = (normStep < tol);
        case 'grad'
            p_k_1=-inv(F.d2f(x_k_1))*F.df(x_k_1);
            alpha=ls(x_k_1,p_k_1,info.alphas(nIter+1));
            x_k=x_k_1+alpha.*p_k_1;
            stopCond = (norm(F.df(x_k), 'inf') < tol*(1 + tol*abs(F.f(x_k))));
    end
    nIter=nIter+1;
    info.xs(:,nIter+1)=x_k;
    info.alphas=[info.alphas alpha]; 
end
nIter=nIter-1;
% Assign output values
xMin = x_k;
fMin = F.f(x_k);

end



