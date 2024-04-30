% show example timewarp procedure
h=[0.01 0.04 0.09 0.16 0.24 0.33 0.42 0.53 0.63 0.72 0.81 0.88 0.94 0.98 1 1 0.98 0.94 0.88 0.81 0.72 0.63 0.53 0.42 0.33 0.24 0.16 0.09 0.04 0.01];
T=convn(randn(100,200),h'*h,'same');
S=T(:,ceil(min(1,(1:.5:size(T,2))/size(T,2)).^.5*size(T,2)));
S=S+convn(.25*randn(size(S)),h'*h,'same');
[F,K] = timewarp(S,T,'display',true);
