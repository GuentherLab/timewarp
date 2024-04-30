function [F,K]=timewarp(S,T,varargin)
% TIMEWARP Finds optimal timewarp between input and output matrices
%
% F = timewarp(S,T)
% finds the function F(t) : [1, length(T)] -> [1, length(S)]
% between target and source time indexes approximating the 
% time-warped source S(:,F(t)) to the target T(:,t)
%
% [F,K] = timewarp(S,T)
% returns also scaling factor K such that: T(:,t) ~= S(:,F(t)) * K(t)
%
% Inputs
%   T [M x Nt] is the target (with M observations and Nt timepoints)
%   S [M x Ns] is the source that we want to align to the target (with M observations and Ns timepoints)
%              (note: if S or T are cell arrays the corresponding cell elements are concatenated along the time (second) dimension [S{:}])
% Outputs
%   F [1 x Nt] is the timewarp vector, with non-decreasing integer values between 1 and Ns
%   K [1 x Nt] is the amplitude scaling vector, with values between .5 and 2 (between 1/cost_capped and cost_capped in general)
%
%
% F = timewarp(S,T, option_name1, option_value1, option_name2, option_value2, ...)
% specifies additional optional arguments. Valid option names are:
%   display         : [true/false, default to false] Displays waitbar and image of timewarp results
%   cost_a          : [0-inf, default 1]            Note: These parameters define cost function parameters. Total cost is computed as:
%   cost_alpha      : [-inf:inf, default 2]             TotalCost = Sum_t{ ...
%   cost_b          : [0-inf, default 0.1]                                 cost_a * distance(Source(F(t))-Target(t))^cost_alpha  + ...
%   cost_beta       : [-inf:inf, default 2]                                cost_b * |(F(t)-F(t-1)) - Ns/Nt|^cost_beta
%   cost_distance   : ['euclidean' 'capped', default 'capped']
%                         For 'euclidean', the distance function is defined as      distance(a,b) = sqrt(mean(abs(a-b).^2))
%                         For 'capped', the distance function is defined as         distance(a,b) = capdist(a,b)
%   cost_capped     : [1-inf, default 2] maximum scaling factor allowed when approximating T~=K*S (with 1/cost_capped <= K <= cost_capped) (only used when cost_distance=='capped')
%   norm_input      : [true/false, default true] Normalizes values in S and T matrices by dividing by the constant term mean(std(T,0,1)) (note: this serves mainly to have the cost_a value scale with the range of values in T)
%
% example:
%   h=[0.01 0.04 0.09 0.16 0.24 0.33 0.42 0.53 0.63 0.72 0.81 0.88 0.94 0.98 1 1 0.98 0.94 0.88 0.81 0.72 0.63 0.53 0.42 0.33 0.24 0.16 0.09 0.04 0.01];
%   T=convn(randn(100,200),h'*h,'same');
%   S=T(:,ceil(min(1,(1:.5:size(T,2))/size(T,2)).^.5*size(T,2)));
%   S=S+convn(.25*randn(size(S)),h'*h,'same');
%   [F,K] = timewarp(S,T,'display',true);

% alfnie@bu.edu
% 4/24
%

if iscell(S), S=cat(2,S{:}); end
if iscell(T), T=cat(2,T{:}); end
Nt=size(T,2);
Ns=size(S,2);
Nd=size(T,1);
assert(size(T,1)==size(S,1),'mismatched number of rows in S and T inputs');
options=struct(...
    'display',false,...
    'norm_input',true,...        % normize S&T to the average <std(T(:,t))> = 1
    'cost_a', 1, ...            % cost = cost_a * |S-T|^cost_alpha   +   cost_b * |dF|^cost_beta
    'cost_alpha', 2, ...
    'cost_b', .1, ...
    'cost_beta', 2,...
    'cost_distance','capped',...
    'cost_capped', 2);

for n=1:2:numel(varargin)
    assert(isfield(options,lower(varargin{n})),'unrecognized option %s',varargin{n});
    options.(lower(varargin{n}))=varargin{n+1};
end
assert(ismember(options.cost_distance,{'euclidean','capped'}),'unrecognized cost_distance value %s',options.cost_distance);
iseuclidean=strcmp(options.cost_distance,'euclidean');

% Forward search
E=nan(Ns,Nt);
E(1,1)=0;
KK=ones(Ns,Nt);
[nil,KK(1,1)]=capdist(T(:,1),S(:,1),options.cost_capped);
IDX=nan(Ns,Nt);
if options.norm_input, norm_input=mean(std(T,0,1)); T=T/norm_input; S=S/norm_input; end
if options.display, hmsg=waitbar(0,'Processing'); end
cost_b_all = options.cost_b*(abs((0:Ns-1)-Ns/Nt).').^options.cost_beta;
%zerosN=zeros(1,N); onesN=ones(1,N);
for nt=2:Nt,
    %E(:,nt)=options.cost_a*sum(abs(S-T(:,nt)).^2,1).';
    for ns=1:Ns,
        % computes source-target difference
        if iseuclidean, 
            cdist=mean(abs(S(:,ns)-T(:,nt)).^2);
        else 
            [cdist,kdist]=capdist(T(:,nt),S(:,ns),options.cost_capped);
            KK(ns,nt)=kdist;
        end
        cost_a = options.cost_a*cdist.^(options.cost_alpha/2);
        cost_b = flipud(cost_b_all(1:ns)); %options.cost_b*abs((1:ns)'-ns).^options.cost_beta;
        [mincost,idxcost]=min(E(1:ns,nt-1)+cost_b);
        E(ns,nt)=mincost+cost_a;
        IDX(ns,nt)=idxcost;
    end
    if options.display, waitbar(nt/Nt,hmsg); end
end
if options.display, close(hmsg); end

% Back-trace the indices list
F=zeros(1,Nt);
F(end)=Ns; %[Energy,F(end)]=min(E(:,end));
for nt=Nt-1:-1:1,
    F(nt)=IDX(F(nt+1),nt+1);
end
K=KK(F+Ns*(0:Nt-1));

if options.display
    idx=findobj('tag',mfilename);
    if isempty(idx), figure('name','timewarp','numbertitle','off','tag',mfilename); else figure(idx); end
    [uE,nill,iE]=unique(E); iE=reshape(iE,size(E));
    clf;
    subplot(221); imagesc(iE); hold on; plot(F,'k.-'); hold off; xlabel('Time (target)'); ylabel('Time (source)'); set(gca,'ydir','normal');
    subplot(223); plot(K,'k.-'); xlabel('Time (target)'); ylabel('Amplitude scaling factor'); if options.cost_capped>1, set(gca,'ylim',[1/options.cost_capped, options.cost_capped]); end
    h1=subplot(322); imagesc(real(S)); title('Source'); set(gca,'xtick',unique(F(10:10:size(T,2))),'xticklabel',[],'fontsize',8); grid on
    h2=subplot(324); imagesc(real(K.*S(:,F))); title('Source time-warped to match Target'); set(gca,'xtick',10:10:size(T,2),'fontsize',8); grid on
    subplot(326); imagesc(real(T)); title('Target'); set(gca,'xtick',10:10:size(T,2),'fontsize',8); grid on; xlabel('Time (target)'); 
    h3=axes('position',[get(h1,'position') get(h2,'position')]*[0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 1; 0 0 1 0 0 0 0 0; 0 1 0 0 0 -1 0 -1]');
    plot([get(h2,'xtick')/diff(get(h2,'xlim'));get(h1,'xtick')/diff(get(h1,'xlim'))],repmat((1:2)',1,numel(get(h1,'xtick'))),'k-','color',.75*[1 1 1]);
    set(gca,'xlim',[0 max([get(h1,'xlim') get(h2,'xlim')])]); 
    axis tight off;
    set(gcf,'color','w');
end


