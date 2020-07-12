function [y] = Rastrigin(x)
% this is a benchmark optimization function
A=10;
[n,runs] = size(x);
m=0;

for i = 1 : n
    m = m + x(i,1)^2 -  A*cos(2*pi*x(i,1));
end

y = 10*n + m;



end