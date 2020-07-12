function model = saeCreate(inputDim, outputDim, options)

% DAECREATE Deep Autoencoder model.
%
%	Description:
%
%	MODEL = DAECREATE(INPUTDIMENSION, OUTPUTDIM, OPTIONS) creates a
%	structure for a multi-layer perceptron. For models with a single
%	hidden layer it is a wrapper structure for NETLAB's multi-layer
%	perceptron model.
%	 Returns:
%	  MODEL - model structure containing the neural network specified.
%	 Arguments:
%	  INPUTDIMENSION - dimension of input data.
%	  OUTPUTDIM - dimension of target data.
%	  OPTIONS - options structure. The structure contains the type of
%	   output 'activation function', the number of hidden units and the
%	   optimiser to be used. A set of default options are given by the
%	   file mlpOptions.
%	
%
%	See also
%	MLPOPTIONS, MLP


%	Copyright (c) 2005, 2006, 2007 Neil D. Lawrence
% 	mlpCreate.m CVS version 1.4
% 	mlpCreate.m SVN version 24
% 	last update 2007-11-03T14:24:25.000000Z

%  Setup and train a stacked denoising autoencoder (SDAE)
rand('state',0)
architecture  = [inputDim options.hiddenDim outputDim];
sae = saesetup(architecture);
sae.ae{1}.activation_function       = options.activation_function;
sae.ae{1}.learningRate              = options.learningRate;
sae.ae{1}.inputZeroMaskedFraction   = options.inputZeroMaskedFraction;
opts.numepochs =   1;
opts.batchsize = size(options.x_tr,2);
sae = saetrain(sae, options.x_tr, opts);
% visualize(sae.ae{1}.W{1}(:,2:end)')

% Use the SDAE to initialize a FFNN
nn = nnsetup(architecture);
nn.activation_function              = options.activation_function;
nn.learningRate                     = options.learningRate;
numParams = 0;
for i = 1:length(architecture)-1
    nn.W{i} = sae.ae{i}.W{1};
    numParams = numParams + numel(nn.W{i});
end

% % Train the FFNN
% opts.numepochs =   1;
% opts.batchsize = 100;
% nn = nntrain(nn, train_x, train_y, opts);
% [er, bad] = nntest(nn, test_x, test_y);

model.type = 'sae';
model.hiddenDim = options.hiddenDim;
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numParams = numParams;
model.outfn = options.output;
model.nn = nn;
model.opts = opts;

model.optimiser = options.optimiser;
