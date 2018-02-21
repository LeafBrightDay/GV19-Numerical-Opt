clear ;
% For computation define as function of 1 vector variable
F.f = @(x) (x(1) - 3*x(2)).^2 + x(1).^4;
F.df = @(x) [2*(x(1) - 3*x(2)) + 4*x(1).^3; -6*(x(1) - 3*x(2))];
F.d2f = @(x) [2 + 12*x(1).^2, -6; -6, 18];
                     
% Parameters
maxIter = 200; 
tol = 1e-10; % Stopping tolerance on relative step length between iterations
debug = 0; % Debugging parameter will switch on step by step visualisation of quadratic model and various step options

% Starting point
x0 = [10; 10]; 

% Trust region parameters 
eta = 0.1;  % Step acceptance relative progress threshold
Delta = 1; % Trust region radius

% Minimisation with 2d subspace and dogleg trust region methods
Fsr1 = rmfield(F,'d2f');
tic;
[xTR_SR1, fTR_SR1, nIterTR_SR1, infoTR_SR1] = trustRegion(Fsr1, x0, @solverCM2dSubspaceExt, Delta, eta, tol, maxIter, debug);
disp(toc/( nIterTR_SR1-1));

err=[];
nIter=[];
for loop=1:(nIterTR_SR1)
    Bk=infoTR_SR1.B{loop};
    xk=infoTR_SR1.xs(:,loop);
    nIter=[nIter loop];
    err(loop)=norm(Bk-F.d2f(xk));
end
plot(nIter,err);
% result = [];
% %plot(infoTR_SR1.xs(1,:), infoTR_SR1.xs(2,:))
% for i = 1:length(infoTR_SR1.xs(1,:))
% result = [result F.f(infoTR_SR1.xs(:,i))];
% end
% plot(1:length(result), result)