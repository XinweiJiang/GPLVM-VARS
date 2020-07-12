function model = saeExpandParam(model, params)

% MLPEXPANDPARAM Update mlp model with new vector of parameters.
%
%	Description:
%
%	MODEL = MLPEXPANDPARAM(MODEL, PARAMS) takes a vector of MLP weights
%	and places them in their respective positions in the MLP model. For
%	single hidden layer neural networks the function is a wrapper for
%	the mlpunpak command.
%	 Returns:
%	  MODEL - the model with the weights distributed in the correct
%	   places.
%	 Arguments:
%	  MODEL - the model in which the weights are to be placed.
%	  PARAMS - a vector of the weights to be placed in the model.
%	
%
%	See also
%	MLPUNPAK, MLPCREATE, MLPEXTRACTPARAM


%	Copyright (c) 2006, 2007 Neil D. Lawrence
% 	mlpExpandParam.m CVS version 1.5
% 	mlpExpandParam.m SVN version 24
% 	last update 2009-09-05T21:46:30.000000Z


startVal = 1;
endVal = (model.inputDim+1)*model.hiddenDim(1);
model.w{1} = reshape(params(startVal:endVal, model.inputDim+1, ...
                        model.hiddenDim(1)));
model.nn.W{1} = model.w{1};

for i = 2:length(model.hiddenDim)
startVal = endVal + 1;
endVal = endVal + (model.hiddenDim(i-1)+1)*model.hiddenDim(i);
model.w{i} = reshape(params(startVal:endVal), model.hiddenDim(i-1)+1, ...
                            model.hiddenDim(i));
model.nn.W{i} = model.w{i};
end
i = length(model.hiddenDim);
startVal = endVal + 1;
endVal = endVal + (model.hiddenDim(i)+1)*model.outputDim;
model.w{i+1} = resphape(params(startVal:endVal), model.hiddenDim(i)+1, ...
                      model.outputDim);
model.nn.W{i+1} = model.w{i+1};

