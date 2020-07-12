
% DEMOILFGPLVM1 Oil data with fully independent training conditional.
%
%	Description:
%	% 	demOilFgplvm1.m SVN version 536
% 	last update 2009-09-29T21:32:46.000000Z

clear all;
clc;
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'USPS';
experimentNo = 1;

% load data
% [Y, lbls] = lvmLoadData(dataSetName);

[x, y, xx, yy] = loadBinaryUSPS(3,5);
xx = [x(51:717,:); xx];
yy = [y(51:717,:); yy];
x1 = x(1:50,:);x2 = x(718:767,:);
x=[x1;x2];
y1 = y(1:50,:);y2 = y(718:767,:);
y=[y1;y2];
clear x1 x2 y1 y2;
xx = x;
yy = y;
Y = x;
lbls = y;


% Set up model
options = fgplvmOptions('ftc');%fitc
options.optimiser = 'optimiMinimize';%scg
latentDim = 2;
d = size(Y, 2);

model = fgplvmCreate(latentDim, d, Y, options);

% Optimise the model.
iters = -100; %1000
display = 1;

model = fgplvmOptimise(model, display, iters);

model.Z = model.X(1:100,:);
model.Y = lbls;
zz = model.X(100:size(model.X,1),:);
model.XX = xx;
model.YY = yy;

zplusY = [model.Z y];
[resultClass, classes, distance] = kNN_SGPLVM(zplusY, zz, 10, model);
result = resultClass - model.YY;
res = tabulate(result)

filename = ['Fgplvm' dataSetName num2str(latentDim)];
plotZ(model, filename);
save(['dem' filename]);
