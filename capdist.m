function [cd,k] = capdist(a,b,lambda)
% d=CAPDIST(a,b)
% mean square distance between a and b pseudo-normalized vectors (normalized up to a scaling factor of 2)
%
% example:
%   x=[1;2]; 
%   K=10;
%   [t1,t2]=ndgrid(linspace(-K,K,1e3), linspace(-K,K,1e3)); 
%   y=[t1(:) t2(:)]'; 
%   d = arrayfun(@(n)capdist(x,y(:,n)),1:size(y,2)); 
%   clf; 
%   contourf(linspace(-K,K,1e3),linspace(-K,K,1e3),reshape(d,[1e3,1e3]),0:1:10);
%   hold on; plot(x(1),x(2),'ko','facecolor','k'); hold off;
%   grid on
%   colormap jet
%   set(gca,'clim',[0 1]);

if nargin<3||isempty(lambda), lambda=2; end

a=a(:);
b=b(:);
na=norm(a);
nb=norm(b);
k=max(1/lambda,min(lambda, na/nb));
sk=sqrt(k);
cd = mean(abs( sk*b-a/sk).^2);
%db=mean(abs(b*k - a).^2);
%da=mean(abs(a/k - b).^2);
%cd = max(da,db);

