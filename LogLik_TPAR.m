function [loglik] = LogLik_TPAR(theta, y_tr, r)

%%% Get Tt
Tt = size(y_tr,1);
dloglik = zeros(Tt,1);

%%% Set Filters
log_lambda_t = zeros(Tt,1);
lambda_t = zeros(Tt,1);

%%% Parameters
omega_1 = theta(1);
beta_1 = theta(2);
alpha_1 = theta(3);
omega_2 = theta(4);
beta_2 = theta(5);
alpha_2 = theta(6);

%%% Initial values
% log_lambda_t(1) = omega/(1-beta);
lambda_t(1) = mean(y_tr,1);

% if y(1) >= r
% log_lambda_t(1) = omega_1/(1-beta_1);
% lambda_t(1) = exp(log_lambda_t(1));
% else
% log_lambda_t(1) = omega_2/(1-beta_2);
% lambda_t(1) = exp(log_lambda_t(1));
% end

if y_tr(1) >= r
log_lambda_t(2) = omega_1/(1-beta_1);
lambda_t(2) = exp(log_lambda_t(2));
else
log_lambda_t(2) = omega_2/(1-beta_2);
lambda_t(2) = exp(log_lambda_t(2));
end

% dloglik(1) = y_t(1)*log(lambda_t(1)) - lambda_t(1);
% loglik = dloglik(1)

dloglik(1:2) = y_tr(1:2).*log(lambda_t(1:2)) - lambda_t(1:2);
loglik = sum(dloglik(1:2));

    %%% Updating
    for t = 2:(Tt-1)
        
        inno_t = (y_tr(t) - lambda_t(t))/(lambda_t(t)); 
        
        if y_tr(t-1)>=r
        log_lambda_t(t+1) = omega_1 + beta_1*log_lambda_t(t) + alpha_1*inno_t;
        else
        log_lambda_t(t+1) = omega_2 + beta_2*log_lambda_t(t) + alpha_2*inno_t;
        end
        
    lambda_t(t+1) = exp(log_lambda_t(t+1));

    dloglik(t+1) = y_tr(t+1)*log(lambda_t(t+1)) - lambda_t(t+1);
    loglik     = loglik + dloglik(t+1); 
    
    end
    loglik = -loglik;
end
