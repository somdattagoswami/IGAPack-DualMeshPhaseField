function [positive_elast]=localDStress2(e1,e2,e3,lambda,miu)
% Computes the positive elastic strain energy and decomposed stress
M = sqrt((e1-e2)^2+4*e3^2);
lambda1 = 0.5*(e1+e2)+0.5*M;
lambda2 = 0.5*(e1+e2)-0.5*M;
positive_elast=(0.5*lambda*(lambda1+lambda2)^2*H(lambda1+lambda2)+miu*(lambda1^2*H(lambda1)+lambda2^2*H(lambda2)));
end

function r=H(a)

if a>=0
    r=1;
else
    r=0;
end
end