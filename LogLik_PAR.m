function [loglik] = LogLik_PAR(theta, y_t)

%%% Get Tt
Tt = size(y_t,1);
dloglik = zeros(Tt,1);

%%% Set Filters
log_lambda_t = zeros(Tt,1);
lambda_t = zeros(Tt,1);

%%% Parameters
omega = theta(1);
beta = theta(2);
alpha = theta(3);

%%% Initial values
% log_lambda_t(1) = omega/(1-beta);
% lambda_t(1) = exp(log_lambda_t(1));
% 
log_lambda_t(1:2) = omega/(1-beta);
lambda_t(1:2) = exp(log_lambda_t(1:2));

% dloglik(1) = y_t(1)*log(lambda_t(1)) - lambda_t(1);
% loglik = dloglik(1)

dloglik(1:2) = y_t(1:2).*log(lambda_t(1:2)) - lambda_t(1:2);
loglik = sum(dloglik(1:2));

    %%% Updating
    for t = 2:(Tt-1)
        
        inno_t = (y_t(t) - lambda_t(t))/(lambda_t(t)); 
        
        log_lambda_t(t+1) = omega + beta*log_lambda_t(t) + alpha*inno_t;
        lambda_t(t+1) = exp(log_lambda_t(t+1));
        
           dloglik(t+1) = y_t(t+1)*log(lambda_t(t+1)) - lambda_t(t+1);
           loglik     = loglik + dloglik(t+1); 
    end
    loglik = -loglik; 
end
