function parameterest

function C=kinetics(theta,t)
m = .0018;k=60.03;b=10000;f=0.01939;e=5000000./2;lambda=20000;lambda2=500;

c0=[491,2000,300000,600000];

[T,Cv]=ode45(@DifEq,t,c0);

function dC=DifEq(t,c)
dcdt=zeros(4,1);
dcdt(1)= theta(3).*k.*c(2).*(c(4))/(c(4)+e)-(theta(5)+m+theta(4))*c(1);
dcdt(2)= theta(1).*k.*c(2).*(c(3))/(c(3)+b)-(f+theta(2)).*c(2);
dcdt(3)= lambda-theta(1).*k.*c(2).*(c(3))/(c(3)+b)-(f).*c(3);
dcdt(4)=lambda2 -theta(3).*k.*c(2).*(c(4))/(c(4)+e)+theta(5).*c(1)-m.*c(4);
dC=dcdt;
end
C=Cv(:,1:1);
end
clear all
clc
close all
t=[0
   1
   2
   3
   4
   5
   6
   7
   8
   9
   
];
c=[491
   220
   302
   375
   712
   927
   773
   344
   319
   354
]; 
theta0=[.011;.001;.001;0.01;.5];
lb=[0.01,.15,0.00001,0.02541,0.35];
ub=[.1,.3,.00031,.423,.7];
[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@kinetics,theta0,t,c,lb,ub);
fprintf(1,'\tRate Constants:\n')
for k1 = 1:length(theta)
fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end
tv = linspace(min(t), max(t)+1);
Cfit = kinetics(theta, tv);

figure(1)
plot(t, c, 'k*','LineWidth',7)
%ylim([0,0.00003]);
hold on
hlp = plot(tv, Cfit,'LineWidth',3);
hold off
xlabel('Time( parameterest')
ylabel('Infected humans')
legend('Actual data','Curve fit', 'Location','N')
end