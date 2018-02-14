function [alpha, info] = backtracking(F, x_k, p, alpha0, opts)
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
% x_k: current iterate
% p: descent direction
% alpha0: initial step length 
% opts: backtracking specific option structure with fields
%   - rho: in (0,1) backtraking step length reduction factor
%   - c1: constant in sufficient decrease condition f(x_k + alpha_k*p) > f(x_k) + c1*alpha_k*(df_k'*p)
%         Typically chosen small, (default 1e-4). 
%
% OUTPUTS
% alpha: step length
% info: structure with information about the backtracking iteration 
%   - alphas: step lengths history 
%
% Copyright (C) 2017  Marta M. Betcke, Kiko Rullan


% ====================== YOUR CODE HERE ======================
info.alphas=[];
while F.f(x_k+alpha0.*p)>F.f(x_k)+opts.c1.*alpha0.*(F.df(x_k))'*p
    alpha0=opts.rho.*alpha0;
    info.alphas=[info.alphas alpha0];
end
% ============================================================
alpha=alpha0;
end
