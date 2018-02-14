function p = solverCM2dSubspace(F, x_k, Delta)
% SOLVERCM2DSUBSPACE Solves quadratic constraint trust region problem via 2d subspace
% p = solverCM2dSubspace(F, x_k, Delta)
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% x_k: current iterate
% Delta: trust region radius
% OUTPUT
% p: step (direction times lenght)
%
% Copyright (C) 2017 Marta M. Betcke, Kiko Rullan 
g=F.df(x_k);
B=F.d2f(x_k);
p_B=-B\g;
if (g(1)*p_B(2)~=g(2)*p_B(1))
    a=B(1,1)*Delta.^2;
    b=B(1,2)*Delta.^2;
    c=B(2,2)*Delta.^2;
    d=g(1,:)*Delta;
    f=g(2,:)*Delta;
    t=roots([-b+d,2*(a-c+f),6*b,2*(-a+c+f),-b-d]);
    p=Delta.*vertcat(2*t'./(1+t'.^2),(1-t'.^2)./(1+t'.^2));
    result=0.5*sum(B*p.*p,1)+g'*p;
    [~,minindex]=min(result);
    p=p(:,minindex);
else
    if g'*B*g<=0
        tol=1;
    else
        tol=min(norm(g).^3/(Delta*g'*B*g),1);
    end
    p=-tol*Delta/norm(g)*g;
end
