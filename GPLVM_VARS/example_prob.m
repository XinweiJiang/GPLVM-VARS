function [y] = example_prob(x)

% this function is same as presented in research paper by Rana et. al

y(1,:) = sin(x(1,:)) +sin(x(2,:)) + exp(x(1,:).*x(2,:)) + 100.*(sqrt(x(1,:)) + sqrt(x(2,:)));

y(2,:) = sin(x(2,:)) +cos(x(3,:)) + exp(x(2,:).*x(3,:)) + (x(2,:).^2) + (x(3,:).^2);




end