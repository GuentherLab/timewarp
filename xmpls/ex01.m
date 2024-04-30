
% show example capdist values
x=[1;2];
K=10;
[tc,tr]=ndgrid(linspace(-K,K,1e3), linspace(-K,K,1e3));
y=[tr(:) tc(:)]';
d = arrayfun(@(n)capdist(x,y(:,n)),1:size(y,2));
clf;
contourf(linspace(-K,K,1e3),linspace(-K,K,1e3),reshape(d,[1e3,1e3]),0:1:50);
hold on; plot(x(1),x(2),'ko','markerfacecolor','k'); hold off;
hold on; plot(x(1)*[1/2 2],x(2)*[1/2 2],'wo','markerfacecolor','w'); hold off;
hold on; text(x(1),x(2)+.5,'SIGNAL','horizontalalignment','center','fontsize',8); hold off;
hold on; text(x(1)*1/2,x(2)*1/2-.5,'minimum-amplitude signal','horizontalalignment','center','fontsize',8); hold off;
hold on; text(x(1)*2,x(2)*2+.5,'maximum-amplitude signal','horizontalalignment','center','fontsize',8); hold off;
grid on
axis equal off
hold on; plot([0 0 nan max(get(gca,'xlim'))*[-1 1]],[max(get(gca,'ylim'))*[-1 1] nan 0 0],'k-'); hold off
