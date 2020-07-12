
function [input_vect_store] = func_proposal_sample_bayes_opt(str,lb,ub)

% The output of this function are the 10 proposal samples obtained via
% bayesian optimization

load forward_sol

tol_bayes = evalin('base','tol_bayes');

lb = transpose(lb);
ub = transpose(ub);
[bb,num_par] = size(lb);
runs = 5000;   % number of latin hyper cube samples for infill sampling
xn = lhsdesign(runs,num_par);
sample_lhs = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb)));
sample_lhs = transpose(sample_lhs);

%% find 10 best performing members from the existing training data
ymp = y';
input_cmg_mat_train = inpt_train;
input_cmg_mat_train = gadd(input_cmg_mat_train,mean_T1n);
input_cmg_mat_train = mapminmax('reverse',input_cmg_mat_train,t1s{1,1});

input_cmg_mat_train = removeconstantrows('reverse',input_cmg_mat_train,t1rows{1,1});

ymp = gadd(ymp,mean_pn);
ymp_real = mapminmax('reverse',ymp,p1s{1,1});
ymp_real_abs = abs(ymp_real);
[ymp_real_abs_norm,y_p1s] = mapminmax(ymp_real_abs,0,1);

% calculate the mean of training means
[uu,vv] = size(ymp_real_abs_norm);

if uu>1
    ymp_mean_store = mean(ymp_real_abs_norm);
else
    ymp_mean_store = (ymp_real_abs_norm);
end
    
num_fit = 10;

% calculate the minimum of training means
for jj = 1 : num_fit
    [min_ymp_val_1,min_ymp_loc_1(1,jj)] = min(ymp_mean_store);
    ymp_mean_store(1,min_ymp_loc_1(1,jj)) = NaN;
    best_inp(:,jj) = input_cmg_mat_train(:,min_ymp_loc_1(1,jj));
    
end

%%  Perform mutation and crossover on the 10 best members

[uinp,vinp]= size(input_cmg_mat_train);

rand_mat = randi([1 num_fit],uinp,num_fit*500);


% Sampling from these best species

for jj = 1 : num_fit*500
    for ii = 1 : uinp
        species_cross(ii,jj) = best_inp(ii,rand_mat(ii,jj));
    end
end
lb_t = transpose(lb);
ub_t = transpose(ub);
species_mut = repmat(best_inp,1,num_fit*100);

[uov,vov] = size(species_mut);


rand_mut = randi([1 uov],int64(uinp/4),vov);
for jj = 1 : vov
    for ii = 1 : int64(uinp/4)
        
        species_mut(rand_mut(ii,jj),jj) = lb_t(rand_mut(ii,jj),1) + rand()*(ub_t(rand_mut(ii,jj),1) - lb_t(rand_mut(ii,jj),1));
    end
end

[next_species] = sample_gen_func(species_cross,species_mut,[]);

%% check if any of the random infill samples or samples from best members overlap with existing training data

[u_lhs,v_lhs] = size(sample_lhs);
[u_next,v_next] = size(next_species);

[input_cmg_norm,inpmap] = mapminmax(input_cmg_mat_train,0,1);
next_species_norm = mapminmax('apply',next_species,inpmap);
sample_lhs_norm = mapminmax('apply',sample_lhs,inpmap);

tol = tol_bayes;   % tolerance of how close infill samples be to the existing training data samples
ct_del = 0;
ct_del_next = 0;
for j = 1 : v_lhs
    
    if j<= v_next
          diff_mat_next{1,j} = gsubtract(input_cmg_norm,next_species_norm(:,j));
    [u_ii,v_ii] = size(diff_mat_next{1,j});
    
    for coln = 1 : v_ii
        norm_diff_mat_next(1,coln) = mean(abs(diff_mat_next{1,j}(:,coln)));
    end
    
   mean_diff_mat_next(:,j) = sum(abs(diff_mat_next{1,j}));
   min_diff_mat_next(1,j) = min(mean_diff_mat_next(:,j));
   
   if min(norm_diff_mat_next) < tol
       ct_del_next = ct_del_next + 1;
       del_mat_next(1,ct_del_next) = j;
   end
    end
        
    
    diff_mat{1,j} = gsubtract(input_cmg_norm,sample_lhs_norm(:,j));
    [u_ii,v_ii] = size(diff_mat{1,j});
    
    for coln = 1 : v_ii
%         norm_diff_mat(1,coln) = norm(diff_mat{1,j}(:,coln));
 norm_diff_mat(1,coln) = mean(abs((diff_mat{1,j}(:,coln))));
    end
    
   mean_diff_mat(:,j) = sum(abs(diff_mat{1,j}));
   min_diff_mat(1,j) = min(mean_diff_mat(:,j));
   
   if min(norm_diff_mat) < tol
       ct_del = ct_del + 1;
       del_mat(1,ct_del) = j;
   end
end

if ct_del_next > 0
del_mat_next = sort(del_mat_next,'descend');
end

if ct_del > 0
del_mat = sort(del_mat,'descend');
end

sample_lhs_old = sample_lhs;
next_species_old = next_species;
 % Remove the infill samples which overlap with the samples in training
 % data
for j = 1 : ct_del
    sample_lhs(:,del_mat(1,j)) = [];
end   

for j = 1 : ct_del_next
    next_species(:,del_mat_next(1,j)) = [];
end  

T_pred_next = next_species;
T_pred = sample_lhs;

sample_vect = (T_pred);
sample_vect_next = (T_pred_next);
% T_pred = sample_vect;
u=1;

for i = 1 : u
    Input_pred_norm_next = removeconstantrows('apply',sample_vect_next, t1rows{i,1});
    Input_pred_norm_next = mapminmax('apply',Input_pred_norm_next, t1s{i,1});
    
    Input_pred_norm = removeconstantrows('apply',sample_vect, t1rows{i,1});
    Input_pred_norm = mapminmax('apply',Input_pred_norm, t1s{i,1});
end

Input_pred_norm = gsubtract(Input_pred_norm,mean_T1n);
Input_pred_norm_next = gsubtract(Input_pred_norm_next,mean_T1n);


[u_t,v_t] = size(T_train);

for jj = 1 : num_dim
    datastruc_ens = datastruc{jj,1};
    ep = 0.001;
    err_msg = 'error';
    stat_chng = 0;
    while strcmp(err_msg,'error') == 1
        try
            [mu_next s2_next] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), Input_pred_norm_next');
        err_msg = 'good';
        catch nosense
            err_msg = 'error';
            hyp_alt = hyp2{jj,1};
            hyp_alt.cov = hyp_alt.cov.*0.5;
            hyp2{jj,1} = hyp_alt;
            stat_chng = stat_chng+1;
        end
    end
        ymp_test_next(jj,:) = mu_next';
    s2_test_next(jj,:) = s2_next';
    
       ep = 0.001;
    err_msg = 'error';
    stat_chng = 0;
    while strcmp(err_msg,'error') == 1
        try
    [mu s2] = gp(hyp2{jj,1}, datastruc_ens.inf , datastruc_ens.mean, datastruc_ens.covfunc, datastruc_ens.likfunc, x, y(:,jj), Input_pred_norm');
     err_msg = 'good';
        catch nosense
            err_msg = 'error';
            hyp_alt = hyp2{jj,1};
            hyp_alt.cov = hyp_alt.cov.*0.5;
            hyp2{jj,1} = hyp_alt;
            stat_chng = stat_chng+1;
        end
    end
    
    ymp_test(jj,:) = mu';
    s2_test(jj,:) = s2';
    
end

%% Since some LMV's can be negative, first unnormalize all LMV's to and take the absolute of all LMV's and normalize it back to [0 1]
ymp_test_next = gadd(ymp_test_next,mean_pn);
ymp_test_real_next = mapminmax('reverse',ymp_test_next,p1s{1,1});
ymp_test_real_abs_next = abs(ymp_test_real_next);
[ymp_real_abs_norm_next] = mapminmax('apply',ymp_test_real_abs_next,y_p1s);

sdev_test_next = s2_test_next.^0.5;
sdev_test_real_next = mapminmax('reverse',sdev_test_next,p1s{1,1});
sdev_test_real_abs_next = abs(sdev_test_real_next);
[sdev_real_abs_norm_next] = mapminmax('apply',sdev_test_real_abs_next,y_p1s);
s2_real_abs_norm_next = sdev_real_abs_norm_next.^2;


ymp_test = gadd(ymp_test,mean_pn);
ymp_test_real = mapminmax('reverse',ymp_test,p1s{1,1});
ymp_test_real_abs = abs(ymp_test_real);
[ymp_real_abs_norm] = mapminmax('apply',ymp_test_real_abs,y_p1s);

sdev_test = s2_test.^0.5;
sdev_test_real = mapminmax('reverse',sdev_test,p1s{1,1});
sdev_test_real_abs = abs(sdev_test_real);
[sdev_real_abs_norm] = mapminmax('apply',sdev_test_real_abs,y_p1s);

%% start the bayesian optimization using EI and PI criteria

s2_real_abs_norm = sdev_real_abs_norm.^2;

[ei_orig_next] = bayesopt_func_ei(hyp2, ymp_real_abs_norm_next', s2_real_abs_norm_next', Input_pred_norm_next', x, ymp_real_abs_norm_next');
[pi_orig_next] = bayesopt_func_pi(hyp2, ymp_real_abs_norm_next', s2_real_abs_norm_next', Input_pred_norm_next', x, ymp_real_abs_norm_next');


[ei_orig] = bayesopt_func_ei(hyp2, ymp_real_abs_norm', s2_real_abs_norm', Input_pred_norm', x, ymp_real_abs_norm');
[pi_orig] = bayesopt_func_pi(hyp2, ymp_real_abs_norm', s2_real_abs_norm', Input_pred_norm', x, ymp_real_abs_norm');


ei_next = ei_orig_next';
pi_next = pi_orig_next';

ei = ei_orig';
pi = pi_orig';

[uu,vv] = size(ei_next);

if uu>1

ei_sum_next = median(ei_next);
pi_sum_next = median(pi_next);

ei_sum = median(ei);
pi_sum = median(pi);

else
    
    ei_sum_next = (ei_next);  % Calculate EI_overall
pi_sum_next = (pi_next);       % Calculate PI overall

ei_sum = (ei);
pi_sum = (pi);
end
    

ei_sum_store_next = ei_sum_next;
pi_sum_store_next = pi_sum_next;

ei_sum_store = ei_sum;
pi_sum_store = pi_sum;

%% Obtain 5 proposal samples that maximizes EI and PI criteria from infill samples of latin hyper cube

tol = 0.1;
for i = 1 : 5
    
    if i <=2
        
        [maxei(1,i), maxei_loc(1,i)] = max(ei_sum);
        ei_sum(1,maxei_loc(1,i)) = NaN;
        input_vect_store(:,i) = T_pred(:,maxei_loc(1,i));
    elseif i>2 && i<=5
        [maxpi(1,i), maxpi_loc(1,i)] = max(pi_sum);
        pi_sum(1,maxpi_loc(1,i)) = NaN;
        loc_diff = gsubtract(maxei_loc,maxpi_loc(1,i));
        while (min(abs(loc_diff)) < tol  )
            pi_sum(1,maxpi_loc(1,i)) = NaN;
            [maxpi(1,i), maxpi_loc(1,i)] = max(pi_sum);
            loc_diff = gsubtract(maxei_loc,maxpi_loc(1,i));
        end
        pi_sum(1,maxpi_loc(1,i)) = NaN;
        input_vect_store(:,i) = T_pred(:,maxpi_loc(1,i));
    end
end

[uu,vv] = size(ei_sum_next);

if vv>=2
    uplim = 4;
elseif vv<2 && vv>=1
    uplim = 0;
else
    uplim = -5;
end

%% Obtain 5 proposal samples that maximizes EI and PI criteria from infill samples of crossover and mutation

for i = 6 : 6+uplim
    
    if i <=7
        
        [maxei_next(1,i), maxei_loc_next(1,i)] = max(ei_sum_next);
        ei_sum_next(1,maxei_loc_next(1,i)) = NaN;
        input_vect_store(:,i) = T_pred_next(:,maxei_loc_next(1,i));
    elseif i>7 && i<=10
        [maxpi_next(1,i), maxpi_loc_next(1,i)] = max(pi_sum_next);
        pi_sum_next(1,maxpi_loc_next(1,i)) = NaN;
        loc_diff_next = gsubtract(maxei_loc_next,maxpi_loc_next(1,i));
        while (min(abs(loc_diff_next)) < tol  )
            pi_sum_next(1,maxpi_loc_next(1,i)) = NaN;
            [maxpi_next(1,i), maxpi_loc_next(1,i)] = max(pi_sum_next);
            loc_diff_next = gsubtract(maxei_loc_next,maxpi_loc_next(1,i));
        end
        pi_sum_next(1,maxpi_loc_next(1,i)) = NaN;
        input_vect_store(:,i) = T_pred_next(:,maxpi_loc_next(1,i));
    end
end

save test_bayesop.mat

end