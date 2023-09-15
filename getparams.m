function [k,z,p] = getparams(gdc,ghf,wc)
    if (gdc > ghf)
        p = wc * sqrt(1 - 2*ghf^2/gdc^2);
    else
        p = wc * sqrt(1 - 2*gdc^2/ghf^2);
        %p = wc * sqrt((ghf^2-2)/(2*gdc^2/ghf^2-ghf^2));
    end
    k = ghf;
    z = p*gdc/ghf;
end

