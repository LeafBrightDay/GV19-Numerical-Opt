clear;
% Initialisation
n = 1e2;
tol = 1e-12;
maxIter = 1e3;
% Initial point
x0 = zeros(n, 1);

% Matrix definition - Use sparse matrices
b = ones(n, 1);
%A1 = diag(1:n);
%A2 = diag([ones(n-1, 1); 100]);
A3 = -diag(ones(n-1, 1), -1) - diag(ones(n-1, 1), 1) + diag(2*ones(n, 1));
% Define xtrue
% xtrue = zeros(n,1);
% xtrue(floor(n/4):floor(n/3)) = 1;
% xtrue(floor(n/3)+1:floor(n/2)) = -2;
% xtrue(floor(n/2)+1:floor(3/4*n)) = 1/2;

%Alternative xtrue


[V, Lambda] = eig(full(A3));
xtrue = V(:,1:n)*randn(n,1);
ev=diag(Lambda)';
ev=ev';%eigvalues
% Identity operator
M = @(y) y;

% Linear preconditioned Conjugate gradient backtracking line search
[xMin1, nIter1, resV1, infoCG1] = conjugateGradient(A3, A3*xtrue, tol, maxIter, M, x0, xtrue);

% Run learner solution.
norm1 = norm(xMin1-xtrue);
perf=zeros(nIter1,1);
for i=1:nIter1
    perf(i)=norm(infoCG1.xs(:,i)-xtrue).^2;
end

hist(ev,length(ev))
figure;semilogy(perf);
title('performance');
xlabel('nIter');
ylabel('log||x-x*||^2');
% Compare.
%assessVariableEqual('norm1', 0, 'AbsoluteTolerance', 1e-10);