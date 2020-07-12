function [perc_inp_mat,num_out_cell,mean_input_mat,cell_sens_resp] = sens_ana_func_gp(lower_thers,upper_thers)

% this function calculates how much change is caused in output values if
% each input parameter is changed from its minimum to maximum value

load forward_sol

input_cmg_saved = input_cmg_mat;
[u_t,v_t] = size(input_cmg_mat);

[u_tv,v_tv] = size(input_cmg_mat);

mean_saved = mean(input_cmg_mat,2);


mean_input_mat = mean(input_cmg_mat,2);
for i = 1: u_t
    max_input_mat(i,1) = max(input_cmg_mat(i,:));
    min_input_mat(i,1) = min(input_cmg_mat(i,:));
end
mean_input_norm = removeconstantrows('apply', mean_input_mat , t1rows{1,1});
mean_input_norm = mapminmax('apply',  mean_input_norm, t1s{1,1});

mean_input_norm = gsubtract(mean_input_norm,mean_T1n);

[u_inp,v_inp] = size(upper_thers);

for i = 1 : u_inp
    lim(i,1) = lower_thers(i,1);
    lim(i,2) = upper_thers(i,1);
end

[lim_norm,lim_ind] = mapminmax(lim,0,1);
ep = 0.001;
err_msg = 'error';
stat_chng = 0;
for jj = 1 : num_dim
    datastruc_ens = datastruc{jj,1};
      while strcmp(err_msg,'error') == 1
    try
    [mu s2] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), mean_input_norm');
     err_msg = 'good';
    catch nosense
        err_msg = 'error';
        hyp_alt = hyp2{jj,1};
        hyp_alt.cov = hyp_alt.cov.*0.5;
        hyp2{jj,1} = hyp_alt;
        stat_chng = stat_chng+1;
    end
    end
    ann_out_mean(jj,:) = mu';
end
ann_out_mean = gadd(ann_out_mean,mean_pn);
ann_out_mean = mapminmax('reverse',ann_out_mean, p1s{1,1});
result_out_mean =  removeconstantrows('reverse',ann_out_mean, p1rows{1,1});


for i = 1 : u_t
    sens_input_cell{i,1} = mean_saved;
    sens_input_cell{i,1}(i,1) = min_input_mat(i,1);

    
    
    sens_input_cell{i,2} = mean_saved;
    sens_input_cell{i,2}(i,1) = max_input_mat(i,1);
    
    sens_input_norm{i,1} = removeconstantrows('apply',sens_input_cell{i,1} , t1rows{1,1});
    sens_input_norm{i,1} = mapminmax('apply',sens_input_norm{i,1}, t1s{1,1});
    sens_input_norm{i,1} = gsubtract(sens_input_norm{i,1},mean_T1n);
    
    sens_input_norm{i,2} = removeconstantrows('apply',sens_input_cell{i,2} , t1rows{1,1});
    sens_input_norm{i,2} = mapminmax('apply',sens_input_norm{i,2}, t1s{1,1});
    sens_input_norm{i,2} = gsubtract(sens_input_norm{i,2},mean_T1n);
    
ep = 0.001;
err_msg = 'error';
stat_chng = 0;
for jj = 1 : num_dim
    
    datastruc_ens = datastruc{jj,1};
    while strcmp(err_msg,'error') == 1
    try
   [mu s2] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), sens_input_norm{i,1}');
    err_msg = 'good';
    catch nosense
        err_msg = 'error';
        hyp_alt = hyp2{jj,1};
        hyp_alt.cov = hyp_alt.cov.*0.5;
        hyp2{jj,1} = hyp_alt;
        stat_chng = stat_chng+1;
    end
    end
    ann_out_1mat(jj,:) = mu';
    
%     stat_chng = 0;
    err_msg = 'error';
     while strcmp(err_msg,'error') == 1
    try
     [mu s2] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), sens_input_norm{i,2}');
    err_msg = 'good';
    catch nosense
        err_msg = 'error';
        hyp_alt = hyp2{jj,1};
        hyp_alt.cov = hyp_alt.cov.*0.5;
        hyp2{jj,1} = hyp_alt;
%         sens_input_norm{i,2}(i,1) = sens_input_norm{i,2}(i,1) - ep;
        stat_chng = stat_chng+1;
    end
    end
     ann_out_2mat(jj,:) = mu';
end
ann_out{i,1} =  ann_out_1mat;
ann_out{i,2} =  ann_out_2mat;
    
ann_out{i,1} = gadd(ann_out{i,1},mean_pn);
ann_out{i,2} = gadd(ann_out{i,2},mean_pn);

        ann_out{i,1} = mapminmax('reverse',ann_out{i,1}, p1s{1,1});
    result_out{i,1} =  removeconstantrows('reverse',ann_out{i,1}, p1rows{1,1});
    
    ann_out{i,2} = mapminmax('reverse',ann_out{i,2}, p1s{1,1});
    result_out{i,2} =  removeconstantrows('reverse',ann_out{i,2}, p1rows{1,1});
    
        perc_inp_mat(i,1) = (sens_input_cell{i,1}(i,1)- mean_input_mat(i,1))/ mean_input_mat(i,1);
    perc_inp_mat(i,2) = (sens_input_cell{i,2}(i,1)- mean_input_mat(i,1))/ mean_input_mat(i,1);
    
    [j_t,m_t] = size(result_out{i,1});
    for j = 1 : j_t
        
        if result_out{i,1}(j,1)>=result_out_mean(j,1)
            num_out_cell{i,1}(j,1) =  abs(result_out{i,1}(j,1)-result_out_mean(j,1));
        else
            num_out_cell{i,1}(j,1) = - abs(result_out{i,1}(j,1)-result_out_mean(j,1));
        end
        
        if result_out{i,2}(j,1)>=result_out_mean(j,1)
            num_out_cell{i,2}(j,1) =  abs(result_out{i,2}(j,1)-result_out_mean(j,1));
        else
            num_out_cell{i,2}(j,1) = - abs(result_out{i,2}(j,1)-result_out_mean(j,1));
        end
        
        tot_perc_out_chng_cell{i,1}(j,1) = num_out_cell{i,2}(j,1) - num_out_cell{i,1}(j,1);
    end
    
    
    
end

for ii = 1 : 2
real_val_perc_inp_mat(:,ii) = perc_inp_mat(:,ii) .* mean_input_mat;
real_val_perc_inp_mat(:,ii) = real_val_perc_inp_mat(:,ii) + mean_input_mat;
end

perc_inp_mat_norm = mapminmax('apply',real_val_perc_inp_mat,lim_ind);


for i = 1 : u_t
    

        tot_perc_inp_chng(i,1) = (perc_inp_mat_norm(i,2) - perc_inp_mat_norm(i,1))./ (lim_norm(i,2) - lim_norm(i,1));

end



for j = 1 : j_t
    for i = 1 : u_t
        
        cell_sens_resp{j,1}(i,1) = i;
        cell_sens_resp{j,1}(i,2) = tot_perc_out_chng_cell{i,1}(j,1)/(tot_perc_inp_chng(i,1)*num_out_cell{i,2}(j,1));
        nan_check = isnan(cell_sens_resp{j,1}(i,2));
        if nan_check == 1
            cell_sens_resp{j,1}(i,2) = 0;
        end
        
    end

    
    cell_sens_resp{j,1}(:,3) = abs(cell_sens_resp{j,1}(:,2))./ max(abs(cell_sens_resp{j,1}(:,2)));
    cell_sens_resp{j,1}(:,4) = cell_sens_resp{j,1}(:,3)./ sum(cell_sens_resp{j,1}(:,3));
    [cell_sens_resp{j,1}(:,5),cell_sens_resp{j,1}(:,6)] = sort(cell_sens_resp{j,1}(:,4),'descend');
end




save sens_ana_result.mat











end