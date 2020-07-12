clc;clear all; close all; st = fclose('all');
randn('seed', 1e7)
rand('seed', 1e7)
Opt.iters = -1000;
dataSetName = 'OilS80';
filename = ['retGpr' dataSetName 'Error' ];

    load('../../../DataSets/Oil/OilS80.mat');
    x = train(:,1:5);
    y = train(:,7);
    xx = test(:,1:5);
    yy = test(:,7);
    
    % Set up the model
    options = gpOptions('ftc');
    options.optimiser = 'optimiMinimize';%scg
%     options.kern = {'rbfard','white'};

    % Use the full Gaussian process model.
    q = size(x, 2);
    d = size(y, 2);
    model = gpCreate(q, d, x, y, options);

    display = 1;
    iters = Opt.iters;

    model = gpOptimise(model, display, iters);

    [mu, varSigma] = gpPosteriorMeanVar(model, xx);

    diffZ = mu - yy;
    retMSE = sqrt(sum(diffZ.*diffZ)/length(yy));
    

save(filename);
