function [hyp2, mu, s2, datastruc,cov_id] = train_gp_mult(x,y,xs)

% this code uses GPML toolbox in order to train a GP model. A list of
% covariance function is used below and uses Empirical Bayes approach in
% order to find the optimum convariance function along with its parameters.

[run,D] = size(x);

likfunc = @likGauss;              % Gaussian likelihood

L = rand(D,1); sf = rand();
ell = rand(); sf = rand();
al = 2;


for mm = 1 : 9
    L = rand(D,1);
    switch mm
        case 1
            cov = {'covSEisoU'};  hypcov = log([ell]);    % isotropic Gaussian
        case 2
            cov = {'covMaternard',3};  hypcov = log([L;sf]);
        case 3
            cov = {'covMaterniso',3};  hypcov = log([ell;sf]);
        case 4
            cov = {@covNNone}; hypcov = log([ell;sf]);           % neural network
        case 5
            cov = {'covRQiso'};  hypcov = log([ell;sf;al]); % ration. quad.
        case 6
            cov1 = {'covMaterniso',3}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 7
            cov1 = {'covMaternard',3}; hypcov1 = log([L;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 8
            cov1 = {@covNNone}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 9
            cov1 = {@covNNone}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covSum',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
            
    end
    
    covcell{mm,1} = cov;
    hypcov_store = hypcov;
    
    for iiv = 1 : 3
        if iiv ==1 
            hyp0.cov = hypcov_store.* abs((rand()*0.1));
        else
        hyp0.cov = hypcov_store.* abs((rand()*iiv) - 0.5);
        end

        mean = [];
        hyp0.mean = [];
        
        hyp0.lik  = -1;
        covfunc = cov;
        
        [hyp2_cell{iiv,mm},~,fval(iiv,mm)] = minimize(hyp0, @gp, -200, @infGaussLik , mean, covfunc, likfunc, x, y);
        
    end
end

[min_fval_row,min_fval_loc] = min(fval);
[min_fval_col,min_fval_loc_col] = min(min_fval_row);

col_num = min_fval_loc_col;

row_num = min_fval_loc(1,col_num);

cov_id = col_num;

hyp2 = hyp2_cell{row_num,col_num};

covfunc = covcell{col_num,1};

datastruc.inf = {'infGaussLik'};
datastruc.mean = mean;
datastruc.covfunc = covfunc;
datastruc.likfunc = likfunc;

ep = 0.001;
err_msg = 'error';
stat_chng = 0;
[ui,vi] = size(xs);
while strcmp(err_msg,'error') == 1
    try
        [mu s2] = gp(hyp2, @infGaussLik , mean, covfunc, likfunc, x, y, xs);
        err_msg = 'good';
    catch nosense
        err_msg = 'error';
        xs_t = xs';
        
       xs_t = xs_t + ones(vi,ui).*ep;
       xs = xs_t';
        stat_chng = stat_chng+1;
    end
end

end