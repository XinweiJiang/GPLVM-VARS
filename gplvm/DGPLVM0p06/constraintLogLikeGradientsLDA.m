function gX = constraintLogLikeGradientsLDA(model)

% CONSTRAINTLOGLIKEGRADIENTSLDA Returns gradients of loglikelihood
%
%	Description:
%	for LDA constraints
%
%	OPTIONS = CONSTRAINTLOGLIKEGRADIENTSLDA(MODEL) Returns loglikelihood
%	for LDAPos constraint
%	 Returns:
%	  OPTIONS - Returns loglikelihood
%	 Arguments:
%	  MODEL - fgplvm model
%	
%
%	See also
%	CONSTRAINTLOGLIKELIHOOD


%	Copyright (c) 2009 Carl Henrik Ek
% 	constraintLogLikeGradientsLDA.m SVN version 295
% 	last update 2009-03-08T08:54:13.000000Z

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
      
      temp = -(dSb_dX*model.A) + dSw_dX;
      gX(model.indices{j}(k),model.dim(i)) = ...
          gX(model.indices{j}(k),model.dim(i)) - model.lambda * ...
          trace((2./model.N)*model.S_b_inv*temp);
      N_acum = N_acum+1;
    end
  end
end

return
