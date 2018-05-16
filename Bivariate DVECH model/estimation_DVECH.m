

%%%%%%%%%%%%%%%%% This MATLAB file contains code for the estimation of the 
%%%%%%%%%%%%%%%%% bivariate DVECH(1,1) model
%%%%%%%%%%%%%%%%% Stock returns of Google and IBM are considered

%% 0. Clean Workspace and Command Window

clear all   %clear workspace
clc         %clear command window     

%% 1. Obtain returns of Google and IBM

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

%obtain the sample size

T=length(r_gog);

%% 2. plot prices and log-returns

subplot(2,2,1)       % add first plot in 2 x 2 grid
plot(p_gog,'k')
grid minor
xlim([0 length(p_gog)])
title('Google prices')


subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(r_gog,'k')
grid minor
xlim([0 length(r_gog)])
title('Google log-returns')


subplot(2,2,3)       % add third plot in 2 x 2 grid
plot(p_ibm,'k')
grid minor
xlim([0 length(p_ibm)])
title('IBM prices')


subplot(2,2,4)       % add fourth plot in 2 x 2 grid
plot(r_ibm,'k')
grid minor
xlim([0 length(r_ibm)])
title('IBM log-returns')




      
%% 3. Optimization Options

      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-6,... % function value convergence criteria 
                         'TolX',1e-6,... % argument convergence criteria
                         'MaxIter',500); % maximum number of iterations    

%% 4. Initial Parameter Values for the optimization
      
      w11=0.2; % initial values for omega
      w12=0.05;
      w22=0.2;
      
      b11=0.8; % initial values for beta
      b12=0.8;
      b22=0.8;
      
      a11=0.1; % initial values for alpha
      a12=0.05;
      a22=0.1;
      
      % store all parameters in the vector theta_ini

      theta_ini = [w11,w12,w22,b11,b12,b22,a11,a12,a22];
            
%% 5. Parameter Space Bounds
        
      lb=[0.0001,-100,0.001,0,0,0,0,0,0];    % lower bound for theta
      ub=[100,100,100,1,1,1,1,1,1];   % upper bound for theta
      
%% 6. Optimize Log Likelihood Criterion
      
      % fmincon input:
      % (1) negative average log-likelihood function: - llik_fun_DVECH()
      % (2) initial parameter: theta_ini
      % (3) parameter space bounds: lb & ub
      % (4) optimization setup: options
      %  Note: a number of parameter restriction are left empty with []

      % fmincon output:
      % (1) parameter estimates: theta_hat
      % (2) log likelihood function value at theta_hat: ls_val
      % (3) exit flag indicating (no) convergence: exitflag
      
      [theta_hat,llik_val,exitflag]=...
          fmincon(@(theta) - llik_fun_DVECH(theta,x),theta_ini,[],[],[],[],lb,ub,[],options);
      
%% 7. Print Output estimation

display('parameter estimates:')
theta_hat

display('log likelihood value:')
-llik_val*T

display('exit flag:')
exitflag

%% 8. Obtain the estimated conditional covariance matrix

% create a Tx3 matrix that will store the conditional ccovariance matrix

VECHt=zeros(T,3);

% set the initial value of the conditional variance equal to the sample
% covariance

C=cov(x);
VECHt(1,:)=[C(1,1),C(1,2),C(2,2)];

for t = 2:T
    
    VECHt(t,1) = theta_hat(1) + theta_hat(4) * VECHt(t-1,1) + theta_hat(7) * x(t-1,1)^2;
    VECHt(t,2) = theta_hat(2) + theta_hat(5) * VECHt(t-1,2) + theta_hat(8) * x(t-1,1)*x(t-1,2);
    VECHt(t,3) = theta_hat(3) + theta_hat(6) * VECHt(t-1,3) + theta_hat(9) * x(t-1,2)^2;

end



%% 8. Conditional variance matrix plots

% we obtain the conditional standard deviations for IBM and Apple 
% log-returns and their conditional correlation

sd1t = VECHt(:,1).^0.5;
sd2t = VECHt(:,3).^0.5;
corr_t = VECHt(:,2)./(sd1t.*sd2t);

% plot the results

subplot(3,1,1)       % add first plot in 3 x 1 grid
plot(sd1t,'k')
grid minor
axis([0 inf min(sd1t) max(sd1t)])
title('conditional sd Google')

subplot(3,1,2)       % add second plot in 3 x 1 grid
plot(sd2t,'k')
grid minor
axis([0 length(sd2t) min(sd2t) max(sd2t)])
title('conditional sd IBM')

subplot(3,1,3)       % add third plot in 3 x 1 grid
plot(corr_t,'k')
grid minor
axis([0 length(corr_t) min(corr_t) max(corr_t)])
title('conditional correlation')
