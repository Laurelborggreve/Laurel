%%%%%%%%%%%%%%%%% The following code shows how to obtain optimal portfolio
%%%%%%%%%%%%%%%%% weights. In this case we use a CCC model to estimate the
%%%%%%%%%%%%%%%%% conditional covariance matrix

%% 00. Clean Workspace and Command Window

clear all   %clear workspace
clc         %clear command window     

%% 0. Obtain monthly returns of Miscrosft and IBM

% Download Google stock prices from yahoo
msft = fetch(yahoo,'MSFT','Adj Close','01/01/1990','01/01/2017','m');
% Download IBM stock prices from yahoo
ibm = fetch(yahoo,'IBM','Adj Close','01/01/1990','01/01/2017','m');

% Obtain Google stock prices  
p_msft = flipud(msft(:,2));
% Obtain IBM stock prices  
p_ibm = flipud(ibm(:,2));


% obtain log-returns for Miscrosft
r_msft = diff(log(p_msft))*100;
% Obtain log-returns for IBM
r_ibm = diff(log(p_ibm))*100;

% Combine the two series of returns in a single matrix


% Obtain the expected log returns using the sample mean
mu1 = mean(r_msft);
mu2 = mean(r_ibm);

x=[(r_msft-mu1),(r_ibm-mu2)];


%% 1. estimate univariate GARCH for each of the series


%% 1a. Optimization Options

      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-12,... % function value convergence criteria 
                         'TolX',1e-12,... % argument convergence criteria
                         'MaxIter',500); % maximum number of iterations    

%% 1b. Initial Parameter values for the optimization
      

      omega_ini = 0.1;% initial value for omega
      alpha_ini = 0.2;  % initial value for alpha
      beta_ini = 0.7;   % initial value for beta
      
      theta_ini = [omega_ini,beta_ini,alpha_ini];
      
%% 1c. Parameter Space Bounds
        
      lb=[0.00000001,0,0];    % lower bound for theta
      ub=[100,1,1];   % upper bound for theta
      
%% 1d. Optimize univariate log-Likelihood functions
   
% Estimate a univariate GARCH(1,1) for Miscrosft log-returns x(:,1)

     [theta_hat1,llik_val1,exitflag1]=...
          fmincon(@(theta) - llik_fun_GARCH(x(:,1),theta),theta_ini,[],[],[],[],lb,ub,[],options);
      
% Estimate a univariate GARCH(1,1) for IBM log-returns x(:,2)

    [theta_hat2,llik_val2,exitflag2]=...
          fmincon(@(theta) - llik_fun_GARCH(x(:,2),theta),theta_ini,[],[],[],[],lb,ub,[],options);
      
      
      
%% 2. estimation of the correlation matrix R from the residuals

%% 2a. Obtain estimated variances

T = length(x(:,1));

s1=zeros(T,1);
s2=zeros(T,1);

s1(1) = var(x(:,1));
s2(1) = var(x(:,2));

% we run the univariate GARCH recursions evaluated at estimated parameter theta_hat

for t=1:(T-1)
    
    s1(t+1) = theta_hat1(1) + theta_hat1(2)*x(t,1)^2 + theta_hat1(3)*s1(t);

    s2(t+1) = theta_hat2(1) + theta_hat2(2)*x(t,2)^2 + theta_hat2(3)*s2(t);
    
end

%% 2b. Obtain standardized series

e1=x(:,1)./sqrt(s1);
e2=x(:,2)./sqrt(s2);

%% 2c. Estimate the correlation of the residuals

R=corrcoef(e1,e2);

%% 2d. Obtain the conditional covariance between Microsoft and IBM

s12=R(1,2)*sqrt(s1).*sqrt(s2); 

%% 3. Optimal portfolio weights

kt = zeros(T,2); % define matrix that will contain the weights


for t=1:T
    
    % set the conditional mean of the log-returns equal to the sample mean.
    % Note that mu1 is the sample mean of Microsft log-returns instead mu2
    % is the sample mean of IBM log returns.
    %(go back to "0." to see how mu1 and mu2 are defined)
    
	mut=[mu1,mu2];
    
    % Obtain the conditional covariance matrix. Note that the vector "s1" 
    % contains the conditional variance of Microsft, "s2" the conditional 
    % variance of IBM and "s12" the conditional covariance
    
 	SIGMAt=[s1(t),s12(t);s12(t),s2(t)];
    
    % Define the portfolio setting the conditional mean vector "mut" and
    % the conditional covariance matrix "SIGMAt"

    p = Portfolio('AssetMean',mut, 'AssetCovar', SIGMAt);
    p = setDefaultConstraints(p);
    
    % Obtain the optimal weights at time t for the portfolio p
    kt(t,:) = estimateMaxSharpeRatio(p);

end

%% plot the results


subplot(2,1,1)
plot(kt(:,1),'k')
xlim([0 inf])
ylim([0 1])
grid minor
title('Portfolio weight to Microsoft k_{1t}')

subplot(2,1,2)
plot(kt(:,2),'k')
xlim([0 inf])
ylim([0 1])
grid minor
title('Portfolio weight to IBM k_{2t}')

      