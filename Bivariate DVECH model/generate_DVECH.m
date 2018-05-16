%% The following code can be used to simulate from a bivariate DVECH model with p=q=1

%% 1. Set sample size and parameter values

T = 1000;

w11=0.1;
w22=0.2;
w12=0.02;

b11=0.7;
b22=0.7;
b12=0.7;

a11=0.2;
a22=0.2;
a12=0.15;


%% 2. Set initial conditions for the conditional variance matrix

% define a Tx2 matrix "x" that will contain the generated series and a Tx3
% matrix "VECHt" that will contain the conditional covariance matrix in
% VECH form. T is the time series length.

x=zeros(T,2);
VECHt=zeros(T,3);

% we choose the unconditional covariance matrix as initialization

VECHt(1,1)=w11/(1-b11-a11);
VECHt(1,2)=w12/(1-b12-a12);
VECHt(1,3)=w22/(1-b22-a22);

%% 3. Genrate from a DVECH model 
%the following loop generates recursively from the DVECH model

for t = 1:T
    
    SIGMAt=[VECHt(t,1),VECHt(t,2);VECHt(t,2),VECHt(t,3)];

    x(t,:)=mvnrnd(zeros(2,1),SIGMAt);    % generate from the observation equation

    VECHt(t+1,1) = w11 + b11 * VECHt(t,1) + a11 * x(t,1)^2;  % updating equation variance series 1
    VECHt(t+1,2) = w12 + b12 * VECHt(t,2) + a12 * x(t,1)*x(t,2); % updating equation covariance between series 1 and series 2
    VECHt(t+1,3) = w22 + b22 * VECHt(t,3) + a22 * x(t,2)^2; % updating equation variance series 2
    
end
    
    
%% 4. Plot series


subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(x(:,1),'k')     % series 1
axis([0 inf -inf inf])
grid minor
title('Series y_{1t}')


subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(x(:,2),'k')     % series 2
axis([0 inf -inf inf])
grid minor
title('Series y_{2t}')



%% 5. Plot conditional covariance

subplot(2,2,1)       % add first plot in 2 x 2 grid
plot(VECHt(:,1),'k')  % varaince series 1
axis([0 inf 0 inf])
grid minor
title('Conditional variance \sigma^2_{1t}')


subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(VECHt(:,3),'k')  % varaince series 2
axis([0 inf 0 inf])
grid minor
title('Conditional variance \sigma^2_{2t}')

subplot(2,2,3)       % add third plot in 2 x 2 grid
plot(VECHt(:,2),'k')  % covariance between series 1 and series 2
axis([0 inf -inf inf])
grid minor
title('Conditional covariance \sigma_{12t}')

subplot(2,2,4)       % add fourth plot in 2 x 2 grid
plot(VECHt(:,2)./(sqrt(VECHt(:,1)).*sqrt(VECHt(:,3))),'k') % correlation between series 1 and series 2
axis([0 inf -inf inf])
grid minor
title('Conditional correlation \rho_{12t}')

