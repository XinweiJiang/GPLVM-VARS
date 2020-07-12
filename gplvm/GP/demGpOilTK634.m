clc;clear all; close all; st = fclose('all');
randn('seed', 1e7)
rand('seed', 1e7)
Opt.iters = -200;
dataSetName = 'TK634';
filename = ['retGpr' dataSetName 'Error' ];

    load('../../../DataSets/Oil/TK634_new.mat');
    x = train(:,2:7);
    y = train(:,1);
    xx = test(:,2:7);
    yy = test(:,1);
    
    % Set up the model
    options = gpOptions('ftc');
    options.optimiser = 'optimiMinimize';%scg
    options.kern = {'rbf','white','bias'};

    % Use the full Gaussian process model.
    q = size(x, 2);
    d = size(y, 2);
    model = gpCreate(q, d, x, y, options);

    display = 1;
    iters = Opt.iters;

    model = gpOptimise(model, display, iters);

    [mu, varSigma] = gpPosteriorMeanVar(model, xx);

    diffZ = mu - yy;
    retMSE = sqrt(sum(diffZ.*diffZ)/length(yy))
    plot([1:length(mu)],mu);hold on;plot([1:length(yy)],yy);
    

save(filename);
