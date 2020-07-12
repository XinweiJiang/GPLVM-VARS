function options = saeOptions(x_tr)

% MLPOPTIONS Options for the multi-layered perceptron.
%
%	Description:
%
%	OPTIONS = MLPOPTIONS returns the default options for a multi-layer
%	perceptron.
%	 Returns:
%	  OPTIONS - default options structure for Multi-layer peceptron.
%
%	OPTIONS = MLPOPTIONS(NUMHIDDEN)
%	 Returns:
%	  OPTIONS - default options structure for Multi-layer peceptron with
%	   the specified number of hidden units.
%	 Arguments:
%	  NUMHIDDEN - number of hidden units.
%	
%
%	See also
%	MLPCREATE, MLP


%	Copyright (c) 2006 Neil D. Lawrence
% 	mlpOptions.m CVS version 1.4
% 	mlpOptions.m SVN version 24
% 	last update 2007-11-03T14:24:25.000000Z

options.x_tr                        = x_tr;
options.hiddenDim                   = [100];
options.activation_function         = 'sigm';
options.learningRate                = 1;
options.inputZeroMaskedFraction     = 0.5;
options.output                      = 'softmax';
options.optimiser = optimiDefaultOptimiser;