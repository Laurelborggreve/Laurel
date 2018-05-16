function llik = llik_fun_DVECH(theta,x)

w11=theta(1);
w12=theta(2);
w22=theta(3);

b11=theta(4);
b12=theta(5);
b22=theta(6);

a11=theta(7);
a12=theta(8);
a22=theta(9);

dim=size(x);
T=dim(1);

VECHt=zeros(T,3);
llik=0;

C=cov(x);
VECHt(1,:)=[C(1,1),C(1,2),C(2,2)];

for t = 2:T
    
    VECHt(t,1)=w11+b11*VECHt(t-1,1)+a11*x(t-1,1)^2;
    VECHt(t,3)=w22+b22*VECHt(t-1,3)+a22*x(t-1,2)^2;
    VECHt(t,2)=w12+b12*VECHt(t-1,2)+a12*x(t-1,1)*x(t-1,2);

    SIGMAt=[VECHt(t,1),VECHt(t,2);VECHt(t,2),VECHt(t,3)];
    
    llik=llik-0.5*(log(det(SIGMAt))+(x(t,:))*inv(SIGMAt)*(x(t,:).'))/T;
    
end

end
