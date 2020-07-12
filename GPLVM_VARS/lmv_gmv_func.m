function [lmv, gmv] = lmv_gmv_func(out, tar)
% this function calculates LMV and GMV for any set of model outputs.
[num_out,Npt] = size(tar);

    term1 = gsubtract(out,tar);
    tar(tar==0) = 1;
    term1 = gdivide(term1,tar);
    term1 = term1.^2;
    term1 = term1./ Npt;
    term2 = (gsubtract(out,tar))./abs(gsubtract(out,tar));
    if Npt > 1
    term1 = sum(term1);
    term2 = sum(term2);
    end
    term1 = term1.^0.5;

lmv = term1.*term2;
if num_out > 1
gmv = sum(abs(lmv))/num_out;
else
    gmv = abs(lmv);
end



end