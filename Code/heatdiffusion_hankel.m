clear; close all;
load('beam.mat')
m = 0:50;
A = diag(-m.^2 * pi^2);
B = (-1).^(m)';
C = [1, 2 * ones(1,size(m,2)-1)];

sys = ss(A,B,C,0);
% Impulse response
[IR,T]=impulse(sys);   size(T) %= 28681 % Associated Hankel matrix

hank=hankel(IR(1:1000,1),IR(1000:2000,1));
size(hank) % = 1000   1001

% Raw discrete-time system
EH=hank(:,1:1000);AH=hank(:,2:1001);BH=hank(:,1);CH=hank(1,1:1000);

% Singular value decompositions
[u1,s1,v1]=svd([EH AH]);
[u2,s2,v2]=svd([EH; AH]);

n = 1:1000;
figure(1); hold on;
set(gca, 'YScale', 'log', 'XScale', 'log');
p1 = loglog(n, diag(s1)/s1(1,1), 'k-');
p2 = loglog(n, diag(s2)/s2(1,1), 'r--');
grid on
xlabel("r")
ylabel("Singular value, \sigma_r")
legend('s1', 's2')

% % Truncation of raw system
k=17; x=u1(:,1:k); y=v2(:,1:k);
Ek=x'*EH*y; Ak=x'*AH*y; Bk=x'*BH; Ck=CH*y; Bcont=Ek\Bk;
dt=T(2,1); % = 9.332e-05
Acont=logm(Ek\Ak)/dt;
sysk=ss(Acont,Bcont,Ck,0);
poles_reg = pole(sys);
poles_trunc = pole(sysk);
figure;hold; 
plot(real(poles_reg), imag(poles_reg),'k.');
plot(real(poles_trunc), imag(poles_trunc),'rx');
xlabel('real(z)')
ylabel('imag(z)')
legend('sys', 'sysk')


[IR,T]=impulse(sys);
[IRk,T]=impulse(sysk,T);
IRk = real(IRk);
figure; hold
% Plot impulse responses and errors
plot(T, IR, 'k-')
plot(T, IRk, 'r--')
xlabel('t')
ylabel('Amplitude')
legend('h', 'h_k')
title('Hankel Impulse Response')


figure; hold on;
bodeplot(sys, 'k-')
bodeplot(sysk, 'r--')
legend('sys', 'sysk')
hold off;


