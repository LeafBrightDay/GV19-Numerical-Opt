function p = solverCM2dSubspaceExt(F, x_k, Delta)
% SOLVERCM2DSUBSPACEEXT Solves quadratic constraint trust region
% p = solverCM2dSubspace(F, x_k, Delta)
% INPUTS
% F: structure with fields
% - f: function handler
% - df: gradient handler
% - d2f: Hessian handler
% x_k: current iterate
% Delta: trust region radius
% OUTPUT
% p: step (direction times lenght)
%
% Copyright (C) 2017 Marta M. Betcke, Kiko Rullan
% Compute gradient and Hessian
g = F.df(x_k);
B = F.d2f(x_k);
% Eigenvalues of Hessian. If B is large, Lanczos - eigs - should be
lambdasB = eig(B);
lambdaB1 = min(lambdasB); %smallest eigenvalue
% Special cases if B has
% 1) negative eigenvalues
% 2) zero eigenvalues
if min(abs(lambdasB)) < eps % zero eigenvalue(s)
    % Take Cauchy point step
    gTBg = g'*(B*g);
    if gTBg <= 0
        tau = 1;
    else
        tau = min(norm(g)^3/(Delta*gTBg), 1);
    end
    p = -tau*Delta/norm(g)*g;
    return;
elseif lambdaB1 < 0 % negative eigenvalue(s)
    alpha = -1.5*lambdaB1; % shift ensuring that B + alpha*I is p.d.
    B = (B + alpha*eye(length(x_k)));
    pNewt = B\g; %2nd order direction
    if norm(pNewt) <= Delta
        npNewt = pNewt/norm(pNewt);
        v = randn(size(x_k));
        v = v/norm(v);
        v = -0.1*npNewt + 0.1*(v - npNewt*npNewt'*v); % v: v'*pNewt <= 0
        p = -pNewt + v; % ensure ||p|| >= ||pNewt||
        %p = -1.1*pNewt; % ensure ||p|| >= ||pNewt||
        %p = -npNewt*Delta; % ensure ||p|| >= ||pNewt||
        return;
    else
        % Orthonormalize the 2D projection subspace
        V = orth([g, pNewt]);
    end
else %positive eigenvalues
    % Orthonormalize the 2D projection subspace
    V = orth([g, B\g]);
end
% Check if gradient and Newton steps are collinear. If so return
if size(V,2) == 1
    % Calculate Cauchy point
    gTBg = g'*(B*g);
    if gTBg <= 0
        tau = 1;
    else
        tau = min(norm(g)^3/(Delta*gTBg), 1);
    end
    p = -tau*Delta/norm(g)*g;
    return;
end
% Project on V
Bv = V'*(B*V);
gv = V'*g;
a = -Bv\gv;
if a'*a < Delta^2
    % Compute the solution p
    p = V*a;
    return;
end
% Case (2)
[Q, Lambdas] = eig(Bv);
lambdas = diag(Lambdas);
Qg = Q'*gv;
r(5) = Delta^2*lambdas(1)^2*lambdas(2)^2 - Qg(1)^2*lambdas(2)^2 -Qg(2)^2*lambdas(1)^2;
r(4) = 2*Delta^2*lambdas(1)^2*lambdas(2) + 2*Delta^2*lambdas(1)*lambdas(2)^2 -2*Qg(1)^2*lambdas(2) - 2*Qg(2)^2*lambdas(1);
r(3) = Delta^2*lambdas(1)^2 + 4*Delta^2*lambdas(1)*lambdas(2) +Delta^2*lambdas(2)^2 -Qg(1)^2- Qg(2)^2;
r(2) = 2*Delta^2*lambdas(1) + 2*Delta^2*lambdas(2);
r(1) = Delta^2;
% Compute roots of the polynomial and select positive one
rootsR = roots(r);
%rootsR = rootsR(rootsR >= 0);
lambda = min(rootsR(rootsR + min(lambdas) > 0));
%lambda = min(rootsR);
% Compute a from Qa(i) = - 1/(lambdas(i) + lambda) * Qg(i)
a = Q* ( (-1./(lambdas(:) + lambda)) .* Qg);
% Compute the solution p
p = V*a;
% Renormalize to ||p|| = Delta, because the condition number of the
p = Delta/norm(p)*p;