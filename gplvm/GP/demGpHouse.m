clc;clear all; close all; st = fclose('all');
randn('seed', 1e7)
rand('seed', 1e7)
Opt.iters = -100;
dataSetName = 'House';
% filename = ['retGp' dataSetName 'Error' ];


% load data
nTr = 304;

load('Housing.mat');
indTr = unidrnd(size(x,1),nTr,1);
indTe = setdiff([1:length(y)], indTr)';
xx = x(indTe,:);
yy = y(indTe,:);
x = x(indTr,:);
y = y(indTr,:);

% Set up the model
options = gpOptions('ftc');
options.optimiser = 'optimiMinimize';%scg
options.kern = {'rbfard','white'};
% % options.kern{1} = 'mlp';
% options.kern = kernCreate(x, {'rbf', 'white'});
% options.kern.comp{1}.inverseWidth = 20;
% options.kern.comp{2}.variance = 0.01;

% Scale outputs to variance 1.
% options.scale2var1 = true;

% Use the full Gaussian process model.
q = size(x, 2);
d = size(y, 2);
model = gpCreate(q, d, x, y, options);

display = 1;
iters = Opt.iters;

model = gpOptimise(model, display, iters);

[mu, varSigma] = gpPosteriorMeanVar(model, xx);

diffZ = mu - yy;
retMSE = sum(diffZ.*diffZ)/length(yy)

filename = ['demGp' dataSetName 'Tr' num2str(size(y,1)) 'Te' num2str(size(yy,1))];
save([filename]);
