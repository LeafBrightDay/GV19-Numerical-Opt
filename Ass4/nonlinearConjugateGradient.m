function [xMin, fMin, nIter, info] = nonlinearConjugateGradient(F, ls, type, alpha0, x0, tol, maxIter)
% NONLINEARCONJUGATEGRADIENT Wrapper function executing conjugate gradient with Fletcher Reeves algorithm
% [xMin, fMin, nIter, info] = nonlinearConjugateGradient(F, ls, 'type', alpha0, x0, tol, maxIter) 
%
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% ls: handle to linear search function
% type: beta update type {'FR', 'PR'}
% alpha0: initial step length 
% rho: in (0,1) backtraking step length reduction factor
% c1: constant in sufficient decrease condition f(x_k + alpha_k*p_k) > f_k + c1*alpha_k*(df_k')*p_k)
%     Typically chosen small, (default 1e-4). 
% x0: initial iterate
% tol: stopping condition on relative error norm tolerance 
%      norm(x_Prev - x_k)/norm(x_k) < tol;
% maxIter: maximum number of iterations
%
% OUTPUTS
% xMin, fMin: minimum and value of f at the minimum
% nIter: number of iterations 
% info: structure with information about the iteration 
%   - xs: iterate history 
%   - alphas: step lengths history 
%
% Copyright (C) 2017  Kiko Rullan, Marta M. Betcke 

% Parameters
% Stopping condition {'step', 'grad'}
stopType = 'step';

% Initialization
nIter = 0;
normError = 1;
x_k = x0;
df_k = F.df(x_k);
p_k = -df_k;
info.xs = x0;
info.alphas = alpha0;
info.betas = [];
stopCond = 0;

% Loop until convergence or maximum number of iterations
while (~stopCond && nIter <= maxIter)
    %============================ YOUR CODE HERE =========================================
    x_k_1=info.xs(:,nIter+1);    
    alpha=min(ls(x_k_1,p_k,alpha0));
    x_k=x_k_1+alpha*p_k;
    beta=(F.df(x_k)'*F.df(x_k))/(F.df(x_k_1)'*F.df(x_k_1));
   % beta=(F.df(x_k)'*(F.df(x_k)-F.df(x_k_1)))/(norm(F.df(x_k_1)).^2);
    p_k=-F.df(x_k)+beta*p_k;
    nIter=nIter+1;
    info.xs(:,nIter+1)=x_k;
    info.alphas=[info.alphas alpha];
    info.betas=[info.betas beta];
    if F.df(x_k) == 0
        break;
    end
    %=====================================================================================
    
    switch stopType
      case 'step' 
        % Compute relative step length
        normStep = norm(x_k - x_k_1)/norm(x_k_1);
        stopCond = (normStep < tol);
      case 'grad'
        stopCond = (norm(df_k, 'inf') < tol*(1 + abs(F.f(x_k))));
    end
end

% Assign output values 
xMin = x_k;
fMin = F.f(x_k); 

