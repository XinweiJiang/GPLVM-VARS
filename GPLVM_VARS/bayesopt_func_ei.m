function [ei] = bayesopt_func_ei(hyp2, mu, s2, xall, x, y)

%% This function calculates Expected Improvement (EI) criteria to perfom Bayesian Optimization.
% Inputs 
% hyp2 = hypothesis of GP model. Can be left as [].
% mu = mean prediction of infill samples
% s2 = standard deviation for all mean predictions of infill samples
% x = input samples, can be left as []
% y = outputs of training data

%% start procedure for bayesian opt
% find minimum y so far
[run,num_dim] = size(y);
for i = 1 : num_dim
[miny(1,i),minyloc(1,i)] = min(y(:,i));
end
sdev = sqrt(s2);
[uiioj,viowe] = size(xall);

for j = 1: num_dim
for i = 1 : uiioj
    u(i,j) = (miny(1,j) - mu(i,j))/ sdev(i,j);
    phis_u(i,j) = normpdf(u(i,j),mu(i,j),sdev(i,j));
    caphis_u(i,j) = cdf('Normal',u(i,j),mu(i,j),sdev(i,j));
    term1 = sdev(i,j) * u(i,j) * phis_u(i,j);
    term2 = sdev(i,j) * caphis_u(i,j);
    
    ei(i,j) = term1 + term2;
   
end
end


end