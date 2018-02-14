function [xMin, nIter, resV, info] = conjugateGradient(A, b, tol, maxIter, M, x0, xtrue)
% CONJUGATE GRADIENT Solves the implicitely symmetric preconditioned CG
%function [x, flag, relres, iter, resvec, xs, V] = mcg(A, b, tol, maxit, M, x0)
% INPUTS
% A: symmetric matrix
% b: specifies the vector to solve Ax = b
% tol: tolerance for the residual
% maxIter: maximum number of iterations
% M: linear preconditioner matrix
% x0: initial iterate
%
% OUTPUTS
% xMin: solution of the system
% nIter: number of iterations taken
% resV: vector of residuals
% info: structure with information about the iteration
%   - xs: iterate history
%   Notice that the left preconditioned NE (with M inner product) and
%   the right preconditioned NE (with M^-1 inner product) produce the same algorithm, 
%   hence only one version.
%
% Copyright (C) 2017 Marta M. Betcke, Kiko Rul·lan

% A and M can be either matrices or function handlers
if ~strcmp(class(A), 'function_handle')
  afun = @(x) A*x;
else
  afun = A;
end
if ~strcmp(class(M), 'function_handle')
  mfun = @(x) M\x;
else
  mfun = M;
end

%===================== YOUR CODE HERE =============================
b=afun(xtrue);
r0=afun(x0)-b;
y0=mfun(r0);
nIter=0;
r_k=r0;
y_k=y0;
p_k=-y0;
x_k=x0;
info.xs=[x0];
resV=[r0];
stopCond=false;

while (~stopCond && nIter<=maxIter)
    alpha=(r_k'*y_k)/(p_k'*afun(p_k));
    x_k=x_k+alpha*p_k;
    r_k_1=r_k+ alpha*afun(p_k);    
    y_k_1=mfun(r_k_1);
    beta=(r_k_1'*y_k_1)/(r_k'*y_k);
    p_k=-y_k_1+beta*p_k;
    r_k=r_k_1;
    y_k=y_k_1;
    nIter=nIter+1;
    if (norm(r_k_1)/norm(b))<tol
        stopCond=true;
    end
    info.xs(:,nIter+1)=x_k;
    resV=[resV;r_k];
end
xMin=x_k;
%==================================================================
