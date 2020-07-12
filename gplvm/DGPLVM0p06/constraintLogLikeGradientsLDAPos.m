function gX = constraintLogLikeGradientsLDAPos(model)

% CONSTRAINTLOGLIKEGRADIENTSLDAPOS Returns loglikegradients for LDAPos constraint
%
%	Description:
%
%	OPTIONS = CONSTRAINTLOGLIKEGRADIENTSLDAPOS(MODEL) Returns
%	loglikegradients for LDAPos constraint
%	 Returns:
%	  OPTIONS - Returns loglikegradients
%	 Arguments:
%	  MODEL - fgplvm model
%	
%
%	See also
%	CONSTRAINTLOGLIKELIHOOD


%	Copyright (c) 2009 Carl Henrik Ek
% 	constraintLogLikeGradientsLDAPos.m SVN version 911
% 	last update 2010-08-23T05:54:38.000000Z

gX = zeros(model.N,model.q);

for(i = 1:1:length(model.dim))
  N_acum = 1;
  for(j = 1:1:model.numClass)
    N_i = length(model.indices{j});
    for(k = 1:1:length(model.indices{j}))
      dOx_dX = zeros(size(model.Ox));
      dOx_dX(model.dim(i),N_acum) = 1;
      
      dSw_dX = dOx_dX*model.W*model.Ox';
      dSb_dX = dOx_dX*model.B*model.Ox';
      
      temp = dSw_dX;
      gX(model.indices{j}(k),model.dim(i)) = ...
          gX(model.indices{j}(k),model.dim(i)) - model.lambda * ...
          trace((2./model.N)*temp);
      N_acum = N_acum+1;
    end
  end
end

return
