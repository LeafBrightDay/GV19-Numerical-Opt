function [x_k, f_k, k, info] = trustRegion(F, x0, solverCM, Delta, eta, tol, maxIter)
% TRUSTREGION Trust region iteration
% [x_k, f_k, k, info] = trustRegion(F, x0, solverCM, Delta, eta, tol, maxIter)
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% x_k: current iterate
% solverCM: handle to solver to quadratic constraint trust region problem
% Delta: upper limit on trust region radius
% eta: step acceptance relative progress threshold
% tol: stopping condition on minimal allowed step
%      norm(x_k - x_k_1)/norm(x_k) < tol;
% maxIter: maximum number of iterations
% OUTPUT
% x_k: minimum
% f_k: objective function value at minimum
% k: number of iterations
% info: structure containing iteration history
%   - xs: taken steps
%   - xind: iterations at which steps were taken
%   - stopCond: shows if stopping criterium was satisfied, otherwsise k = maxIter
%   
% Reference: Algorithm 4.1 in Nocedal Wright
%
% Copyright (C) 2017 Marta M. Betcke, Kiko Rullan 
k=0;
x_k=x0;
delta=0.1;
info.xs=x0;
info.xind(k+1)=k;
stopCond=false;

while (~stopCond && k <= maxIter)
    x_k_1=info.xs(:,k+1);
    p=solverCM(F,x_k_1,delta);
    rho=(F.f(x_k_1)-F.f(x_k_1+p))/(-F.df(x_k_1)'*p-0.5*(p')*F.d2f(x_k_1)*p);
    if rho<0.25
        delta=0.25*norm(p);
    else
        if rho>0.75 && (norm(p)==delta)    
            delta=min(2*delta,Delta);
        end
    end
    if rho>eta
        x_k=x_k_1+p;
    end
    stopCond=(norm(x_k - x_k_1)/norm(x_k) < tol);
    k=k+1;
    info.xind(k+1)=k;
    info.xs(:,k+1)=x_k;   
end
f_k=F.f(x_k);
if stopCond
   info.stopCond=stopCond;
end
if ~stopCond
    k=maxIter;
end