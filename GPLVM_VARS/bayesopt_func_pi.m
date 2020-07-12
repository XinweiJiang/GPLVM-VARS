function [pi] = bayesopt_func_pi(hyp2, mu, s2, xall, x, y)

%% This function calculates another version of Expected Improvement (EI) which we call as (PI) criteria to perfom Bayesian Optimization. This is similar to Probability Improvement (PI) criterion.
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
[uuuu,dsf] = size(xall);
for j = 1: num_dim
for i = 1 : uuuu
    u(i,j) = (miny(1,j) - mu(i,j));
     u(i,j) = (miny(1,j));
     pi(i,j) = u(i,j) - (mu(i,j)-2*sdev(i,j));
%     pi(i,1) = cdf('Normal',u(i,1),mu(i,1),sdev(i,1));
  
end
end


end