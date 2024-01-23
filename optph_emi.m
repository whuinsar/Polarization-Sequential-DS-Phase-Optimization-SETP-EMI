function [slcopt,ph] = optph_emi(covMT,rep)
gamma1 = Regularization(abs(covMT),size(covMT,1),rep);
[U,~] = svd(pinv(gamma1).*covMT);
kU = (U(:,end));
kU = sqrt(size(covMT,1))*kU/norm(kU);
% kU = 1./kU;
tph = angle(kU);
kU = exp(1j*(tph -tph(1)));
ph = angle(kU);
slcopt = kU;

end