
% DEMVOWELSFGPLVM1 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints.
%
%	Description:
%	% 	demVowelsFgplvm1.m SVN version 536
% 	last update 2009-09-29T21:25:56.000000Z

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'vowels';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);
Y = Y(1:200,:);lbls = lbls(1:200,:);

% Set up model
options = fgplvmOptions('ftc');%fitc
options.back = 'mlp';
options.backOptions = mlpOptions;
options.numActive = 200;
latentDim = 2;
d = size(Y, 2);

model = fgplvmCreate(latentDim, d, Y, options);

% Optimise the model.
iters = 10;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

% Load the results and display dynamically.
lvmResultsDynamic(model.type, dataSetName, experimentNo, 'vector')
% lvmResultsDynamic(dataSetName, experimentNo, 'vector')

errors = lvmNearestNeighbour(model, lbls);
