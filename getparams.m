function [k,z,p] = getparams(gdc,ghf,wc)
%p = sqrt(wc^2/(2*gdc^2/ghf^2 - 1));
%
p = 2.4*wc * sqrt(1 - 2*ghf^2/gdc^2);
k = ghf;
z = p*gdc/ghf;
end

