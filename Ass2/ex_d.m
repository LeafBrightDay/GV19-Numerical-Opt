clear;close all;
% For computation define as function of 1 vector variable
F.f = @(x) 100.*(x(2) - x(1)^2).^2 + (1 - x(1)).^2; % function handler, 2-dim vector
F.df = @(x) [-400*(x(2) - x(1)^2)*x(1) - 2*(1 - x(1)); 
              200*(x(2) - x(1)^2)];  % gradient handler, 2-dim vector
F.d2f = @(x) [-400*(x(2) - 3*x(1)^2) + 2, -400*x(1); -400*x(1), 200]; % hessian handler, 2-dim vector
% For visualisation proposes define as function of 2 variables
rosenbrock = @(x,y) 100.*(y - x.^2).^2 + (1 - x).^2;


% Initialisation
alpha0 = 1;
maxIter = 1e4;
alpha_max = alpha0;
tol = 1e-6;
%=============================
% Point x0 = [1.2; 1.2]
%=============================
x0 = [1.2; 1.2];
x1 = [-1.2; 1];

% descent line search
opts.c1=1e-4;
%opts.c2=0.9; %for newton
opts.c2=0.9;
lsFun = @(x_k, p_k, alpha0) lineSearch(F, x_k, p_k, alpha_max, opts);


%lsFun = @(x_k, p_k, alpha0) zoomInt(F, alpha_l, alpha_h, c1, c2);
%[xSteep1, fSteep1, nIterSteep1, infoSteep1] = descentLineSearch(F, 'steepest', lsFun, alpha0, x0, tol, maxIter);
[xSteep1, fSteep1, nIterSteep1, infoSteep1] = descentLineSearch(F, 'steepest', lsFun, alpha0, x1, tol, maxIter);

%[xSteep2, fSteep2, nIterSteep2, infoSteep2] = descentLineSearch(F, 'newton', lsFun, alpha0, x0, tol, maxIter);
[xSteep2, fSteep2, nIterSteep2, infoSteep2] = descentLineSearch(F, 'newton', lsFun, alpha0, x1, tol, maxIter);
n_iter1=zeros(1,nIterSteep1+2);
for loop=1:(nIterSteep1+2)
    n_iter1(loop)=loop-1;
end
n_iter2=zeros(1,nIterSteep2+2);
for i=1:(nIterSteep2+2)
    n_iter2(i)=i-1;
end
subplot(1,2,1);loglog(n_iter1,infoSteep1.alphas);
subplot(1,2,2);plot(n_iter2,infoSteep2.alphas);