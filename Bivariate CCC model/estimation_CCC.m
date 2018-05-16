%%%%%%%%%%%%%%%%% The following code shows how to estimate a bivariate CCC
%%%%%%%%%%%%%%%%% model using the equation by equation approach

%% 00. Clean Workspace and Command Window

clear all   %clear workspace
clc         %clear command window     

%% 0. Obtain returns of Google and IBM

% download Google stock prices from yahoo
gog = fetch(yahoo,'GOOG','Adj Close','01/01/2005','01/01/2017','d');
% download IBM stock prices from yahoo
ibm = fetch(yahoo,'IBM','Adj Close','01/01/2005','01/01/2017','d');

% Obtain Google stock prices  
p_gog = flipud(gog(:,2));
% Obtain IBM stock prices  
p_ibm = flipud(ibm(:,2));

% obtain log-returns for Google
r_gog = diff(log(p_gog))*100;
% obtain log-returns for IBM
r_ibm = diff(log(p_ibm))*100;

%combine the two series of returns in a single matrix
x=[r_gog,r_ibm];

%% 1. estimate univariate GARCH for each of the series

%% 1a. Optimization Options

      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-12,... % function value convergence criteria 
                         'TolX',1e-12,... % argument convergence criteria
                         'MaxIter',500); % maximum number of iterations    

%% 1b. Initial Parameter Values for the optimization
      

      omega_ini = 0.1;% initial value for omega
      alpha_ini = 0.2;  % initial value for alpha
      beta_ini = 0.7;   % initial value for beta
      
      theta_ini = [omega_ini,beta_ini,alpha_ini];
      
%% 1c. Parameter Space Bounds
        
      lb=[0.00000001,0,0];    % lower bound for theta
      ub=[100,1,1];   % upper bound for theta
      
%% 1d. Optimize univariate log-Likelihood functions
   
% Estimate a univariate GARCH(1,1) for Google log-returns x(:,1)

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

%% 3. Display estimation results

%conditional correlation:
R(1,2)

% w1, a1, b1 from variance equation 1st series:
theta_hat1

% w2, a2, b2 from variance equation 2nd series:
theta_hat2

%% 4. Plot estimated conditional variances, covariance and correlation


subplot(2,2,1)
plot(s1,'k')
xlim([0 inf])
ylim([-inf inf])
grid minor
title('Google \sigma^2_{1t}')

subplot(2,2,2)
plot(s2,'k')
xlim([0 inf])
ylim([-inf inf])
grid minor
title('IBM \sigma^2_{2t}')
      

subplot(2,2,3)
plot(s12,'k')
xlim([0 inf])
ylim([-inf inf])
grid minor
title('Covariance \rho^2_{12t}')

subplot(2,2,4)
plot(s12./(sqrt(s1).*sqrt(s2)),'k')
xlim([0 inf])
ylim([0 1])
grid minor
title('Correlation \rho^2_{12t}')
