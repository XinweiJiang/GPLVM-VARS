function [params, names] = saeExtractParam(model)

% MLPEXTRACTPARAM Extract weights and biases from an MLP.
%
%	Description:
%
%	[PARAMS, NAMES] = MLPEXTRACTPARAM(MODEL) returns a vector of all the
%	weights and biases from a multi-layer perceptron model. For single
%	hidden layer models the function is a wrapper for the mlppak
%	command.
%	 Returns:
%	  PARAMS - vector of all the weights and biases returned by the
%	   model. The structure is governed by mlppak.
%	  NAMES - optional additional returned cell array of the names of
%	   the parameters.
%	 Arguments:
%	  MODEL - the model from which we wish to extract the weights and
%	   biases.
%	
%
%	See also
%	MLPPAK, MLPCREATE, MLPEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2006, 2007 Neil D. Lawrence
% 	mlpExtractParam.m CVS version 1.6
% 	mlpExtractParam.m SVN version 24
% 	last update 2009-09-05T21:46:28.000000Z


params = zeros(1, model.numParams);
startVal = 1;
endVal = (model.inputDim+1)*model.hiddenDim(1);
params(startVal:endVal) = model.nn.W{1}(:)';

for i = 2:length(model.hiddenDim)
startVal = endVal + 1;
endVal = endVal + (model.hiddenDim(i-1)+1)*model.hiddenDim(i);
params(startVal:endVal) = model.nn.W{i}(:)';
end
i = length(model.hiddenDim);
startVal = endVal + 1;
endVal = endVal + (model.hiddenDim(i)+1)*model.outputDim;
params(startVal:endVal) = model.nn.W{i+1}(:)';

if nargout > 1
counter = 0;
for j = 1:size(model.w{1}, 2)
  for i = 1:size(model.w{1}, 1)
    counter = counter + 1;
    names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
  end
end

for k = 2:length(model.hiddenDim)
  for j = 1:size(model.w{k}, 2)
    for i = 1:size(model.w{k}, 1)
      counter = counter + 1;
      names{counter} = ['Hidden weight layer ' num2str(k-1) '-' ...
                        num2str(k) ', node '  num2str(i) '-' ...
                        num2str(j)];
    end
  end
end

for j = 1:size(model.w{end}, 2)
  for i = 1:size(model.w{end}, 1)
    counter = counter + 1;
    names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
  end
end
end
