clear all; close all; clc;

x=0:0.01:10;
l1 = levy_density(x,1); l2 = levy_density(x,2); l3 = levy_density(x,3); l4 = levy_density(x,4);

% Levy density functions
figure;
plot(x,l1,'b',x,l2,'r',x,l3,'k',x,l4,'m');
h=xlabel('x'); set(h,'Fontsize',16); h=ylabel('p(x|\sigma)'); set(h,'Fontsize',16);
h=legend('\sigma=1','\sigma=2','\sigma=3','\sigma=4'); set(h,'Fontsize',14);

% Comparison Levy-Rayleigh-Poisson
sigma = 1;
%g1 = 1/(sqrt(2*pi)*sigma) .* exp(-x.^2 / (2*sigma) );
r1 = x/(sigma^2) .* exp(-x.^2 /(2*sigma^2));
p1= poisson_density(x,1);

figure;
subplot(1,2,1);
plot(x,r1,'b--',x,p1,'r-.',x,l1,'k');
h=xlabel('x'); set(h,'Fontsize',16); h=ylabel('p(x|\sigma)'); set(h,'Fontsize',16); 
h=legend('Rayleigh','Poisson','Lévy'); set(h,'Fontsize',14);

subplot(1,2,2);
semilogy(x(100:end),r1(100:end),'b--',x(100:end),p1(100:end),'r-.',x(100:end),l1(100:end),'k');
h=xlabel('x'); set(h,'Fontsize',16); h=ylabel('p(x|\sigma)'); set(h,'Fontsize',16);
aa=axis; axis([x(100) x(end) aa(3) aa(4)]);

