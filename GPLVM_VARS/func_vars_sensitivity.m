  function [max_area_final,max_area_ind_final] = func_vars_sensitivity()
% this function provides local and global sensitivity metrics from VARS
% based sensitivity analysis
  % max_area_final = global relative sensitivity index
  % max_area_ind_final = local relative sensitivity index
  
  load forward_sol   % load trained gpfs model in order to perform sensitivity analysis

  
  % Calciulate mean input vector
mean_input_mat = mean(input_cmg_mat,2);
[u_t,v_t] = size(input_cmg_mat);
for i = 1: u_t
    max_input_mat(i,1) = max(input_cmg_mat(i,:));
    min_input_mat(i,1) = min(input_cmg_mat(i,:));
end

% normalize the minimum, mean and maximum input vectors between [0 1]
min_input_norm = removeconstantrows('apply',min_input_mat , t1rows{1,1});
min_input_norm = mapminmax('apply',  min_input_norm, t1s{1,1});
min_input_norm = gsubtract(min_input_norm,mean_T1n);

max_input_norm = removeconstantrows('apply',max_input_mat , t1rows{1,1});
max_input_norm = mapminmax('apply',  max_input_norm, t1s{1,1});
max_input_norm = gsubtract(max_input_norm,mean_T1n);

mean_input_norm = removeconstantrows('apply', mean_input_mat , t1rows{1,1});
mean_input_norm = mapminmax('apply',  mean_input_norm, t1s{1,1});
mean_input_norm = gsubtract(mean_input_norm,mean_T1n);

targ_norm = mapminmax('apply',comb_inp{1,1},p1s{1,1});
transp_targ = transpose(targ_norm);

[u_tt_targ,v_tt_targ] = size(transp_targ);

% create a random matrix for all input parameters with 100 random entries
run_tot = 100;
for imi = 1 : run_tot
    rand_input(:,imi) = min_input_norm + rand().*(max_input_norm - min_input_norm);
    mean_input_all(:,imi) = mean_input_norm;
end

[u_t,vere] = size(mean_input_norm);

% Perform VARS based sensitivity analysis
v_tt = u_t;
 for tt = 1 : v_tt_targ
    ct = 0;
    for jj = 1 : v_tt
        for ii = 1 : v_tt
                        ct = ct+1;
                  sens_input = mean_input_all;
                  sens_input(ii,:) = rand_input(ii,:);  % vary two inputs at a time to account for second order interactions
                  sens_input(jj,:) = rand_input(jj,:);  % vary two inputs at a time to account for second order interactions
                  transp_input = sens_input';
                  
                  ann_out_sens = zeros(v_tt_targ,run_tot);
                  
                  % calculate output of the gpfs model for a specific set
                  % of input vectors
                  
                  ep = 0.001;
                  err_msg = 'error';
                  stat_chng = 0;
                  for tt_ins = tt
                      datastruc_ens = datastruc{tt_ins,1};
                      while strcmp(err_msg,'error') == 1
                          try
                              [mu s2] = gp(hyp2{tt_ins,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,tt_ins), transp_input);
                              err_msg = 'good';
                          catch nosense
                              err_msg = 'error';
                              hyp_alt = hyp2{tt_ins,1};
                              hyp_alt.cov = hyp_alt.cov.*0.5;
                              hyp2{tt_ins,1} = hyp_alt;
                              stat_chng = stat_chng+1;
                          end
                      end
                      ann_out_sens(tt_ins,:) = mu';
                  end
                  ann_out_sens = gadd(ann_out_sens,mean_pn);
                  
                  transp_targ = ann_out_sens';
                  
                  %% Calculate variogram 
                  if ii == jj
       d{ct,tt} = variogram([transp_input(:,jj) transp_input(:,ii)],transp_targ(:,tt),'plotit',false,'nrbins',50,'maxdist',0.8);
                  else
            d{ct,tt} = variogram([transp_input(:,jj) transp_input(:,ii)],transp_targ(:,tt),'plotit',false,'nrbins',50,'maxdist',0.8);
                  end
                  
                  %% Calculate area under the variogram, IVARS metric
                  dind = d{ct,tt};
        var_area_ind = cumtrapz(dind.distance,dind.val);
        var_area{1,tt}(ii,jj) = var_area_ind(end,1);
        var_area{1,tt}(jj,ii) = var_area{1,tt}(ii,jj);
        

        end
    end
    
    %% sum IVARS of individual input parameters for a specific output or LMV
    var_area{1,tt}(isnan(var_area{1,tt}) > 0) = 0;
       var_area_norm{1,tt} = var_area{1,tt};
    max_area = sum(var_area_norm{1,tt});
    max_area_mat(:,tt) = max_area';
 end

 % Calculate global relative sensitivity index
sum_max_area = sum(max_area_mat,2);
    max_area_final_init =  sum_max_area./ max(max(sum_max_area));
     max_area_final_init = - max_area_final_init ;
     
     % Calculate local relative sensitivity index
     max_area_final_individual_init = - gdivide(max_area_mat,max(max_area_mat));
     max_area_final_init = min(max_area_final_individual_init,[],2);

max_area_final = removeconstantrows('reverse',max_area_final_init,t1rows{1,1});
max_area_ind_final = removeconstantrows('reverse',max_area_final_individual_init,t1rows{1,1});
max_area_final(max_area_final > 0) = 0;
max_area_ind_final(max_area_ind_final > 0) = 0;

max_area_final = -max_area_final;   % global relative sensitivity index
max_area_ind_final = -max_area_ind_final;   % local relative sensitivity index

%% plot graphs of local sensitivity
for j = 1 : v_tt_targ
    figure(100+j)
    bar(max_area_ind_final(:,j));
    title(['local sensitivity w.r.t output' num2str(j)])
    xlabel('Parameters')
    ylabel('Relative Sensitivity')
end

%% plot graph of global sensitivity

figure(100+j+1)
 bar(max_area_final(:,1));
 title(['Global relative sensitivity'])
    xlabel('Parameters')
    ylabel('Relative Sensitivity')
    
    save vars_sens_saved.mat
 
 end