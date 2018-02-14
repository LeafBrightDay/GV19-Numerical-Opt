clear all;
close all

x = linspace(-5,5,100);
y = x;
[X,Y] = meshgrid(x,y);

figure, contour(X, Y, f(X,Y),'ShowText','on')
figure, surf(X, Y, f(X,Y),'EdgeColor','none')

%compare gradient & diff



function result=f(x,y)
result=2.*x+4.*y+x.^2-2*(y.^2);
end