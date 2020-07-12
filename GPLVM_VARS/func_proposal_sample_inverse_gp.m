
function [Target_pred_result] = func_proposal_sample_inverse_gp()

% this function gives one proposal sample via inverse GP

load inverse_model.mat

pred_inp = (comb_inp{1,1}(:,2)./max(comb_inp{1,1}(:,2))).*0.000000001;  % give very low value of LMVs
% pred_inp = comb_inp{1,1}(:,2).*0.001;
u=1;
sample_vect = pred_inp;
for i = 1 : u
%   Input_pred_norm = removeconstantrows('apply',sample_vect, p1rows{i,1});
    Input_pred_norm = mapminmax('apply',sample_vect, p1s{i,1});
end

Input_pred_norm = gsubtract(Input_pred_norm,mean_pn);

Input_pred_norm =(Input_pred_norm);
v_t=1;

[u_t,v_t] = size(T_train);

if bUseGPLVM4GPIS == 0    
    for jj = 1 : num_dim
        datastruc_ens = datastruc{jj,1};
        [mu s2] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), Input_pred_norm');
        ymp_test(jj,:) = mu';
        s2_test(jj,:) = s2';  
    end
else
        % In GPLVMIS we obtain the temporary (or proposal) solution
        % directly without loops.
        options = fgplvmOptions('ftc');%fitc
        options.optimiser = 'optimiMinimize';%scg
        n2 = dist2(x, Input_pred_norm');
        [distmin, minindex] = min(n2);
        
        options.initX = [y; y(minindex,:)];
%         options.oriX_NoLastOne = y;
  
        model = fgplvmCreate(num_dim, num_lvm, [x; Input_pred_norm'], options);% 
        model = fgplvmOptimise(model, 1, -200);
        ymp_test = model.X(model.N, :)';
end


ymp_test = gadd(ymp_test,mean_T1n);

    linreg_pred_result = ymp_test ;
% perform back 
lin_pred_real_1 = mapminmax('reverse',linreg_pred_result, t1s{1,1});
lin_pred_real = removeconstantrows('reverse',lin_pred_real_1, t1rows{1,1});


Target_pred_result = (lin_pred_real);
end