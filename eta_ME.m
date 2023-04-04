clear all; close all; clc

A=10.^; La=length(A);
eta_matlab=zeros(1,La);
eta_newton=zeros(1,La);
eta_dicho=zeros(1,La); eta_DL=zeros(1,La);
eta0=1; niter = 1;

cost=@(eta,a) exp(2.*eta.*log(a))-1-2.*eta.*a.*log(a);

for ia=1:La
    a=A(ia);
    ff = @(eta) cost(eta,a);
    
    % Matlab estimate
     eta_matlab(ia) = fzero(ff,1);
    
    
    % Approx with Newton method
    eta_a = eta0;
    for m=1:niter
        eta_a = eta_a+(1+2*a*log(a)*eta_a-a^(2*eta_a))/(2*a*log(a)*(a^(2*eta_a-1)-1));
    end
    eta_newton(ia)=eta_a;
    
     % Approx with Dichotomy
    [auxm,auxp] = find_zero_dicho(a,niter);
    eta_dicho(ia)=auxm;
    
    % Approx with order-1 DL
   % eta_DL(ia) = 1/2+(a*log(a)-a+1)/(a-1)^2;
   
   eta_a_inf_1(ia) = 1/2+sqrt((a*log(a)-a+1)/(2*a*log(a)^2)); %DL en 1/2
   eta_a_sup_1(ia)= 1/2 + (a.*log(a)-a+1)./(2*(a.^2-a-a.*log(a))); % corde 1/2->1 
   
end
eta_DL = eta_a_inf_1; eta_DL(A>1) = eta_a_sup_1(A>1);

% Plot estimators
figure;
subplot(1,2,1);
plot(A,eta_matlab,'k',A,eta_newton,'b-.',A,eta_DL,'r--'); ha=legend('$\tilde{\eta}$','$\eta^*_{newton}$','$\eta^*_{DL}$');set(ha,'FontSize',14,'interpreter','latex');
ha=xlabel('$a$'); set(ha,'FontSize',16,'interpreter','latex');
ha=ylabel('$\eta$'); set(ha,'FontSize',16,'interpreter','latex');

subplot(1,2,2);
plot(A,100*abs(eta_matlab-eta_newton)./eta_matlab,'b-.',A,100*abs(eta_matlab-eta_DL)./eta_matlab,'r--');
ha=xlabel('$a$'); set(ha,'FontSize',16,'interpreter','latex');
ha=ylabel('Erreur (%)'); set(ha,'FontSize',16);
