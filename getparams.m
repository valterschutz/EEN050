function [k,z,p] = getparams(gdc,ghf,wc)
p = sqrt(wc^2/(2*gdc^2/ghf - 1));
k = ghf;
z = p*gdc/ghf;
end

