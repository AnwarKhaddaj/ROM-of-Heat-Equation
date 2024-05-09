%  Reduced Models of Heat Diffusion Equation
%  ELEC519 Final Project
%
%  AUTHORS:  Karthik Goli, Anwar Khaddaj, and Carlos Taveras
%            Rice University
%            April 12, 2024
%            Last Modified by Anwar: April 16, 2024 3:20pm

% Close all previous figures and clear all variables
close all;
clear all;

%% Plotting impulse responses z(t)
t = linspace(0,1,1000); % time step = 0.1 on [0,1]
z = 1+2*h(t);
z_tilde = 1+2*h_tilde(t);
z_mod = 1+2*(-exp(-pi^2*t)+exp(-4*pi^2*t));

% balanced truncation
n = 0:50;
A = diag(-n.^2 * pi^2);
B = (-1).^(n)';
C = [1, 2 * ones(1,size(n,2)-1)];

sys = ss(A,B,C,0);
%Plotting poles 
poles = pole(sys);
figure;
plot(real(poles),imag(poles),'.');
title('Poles of the original system');

[sysb,sv,T]=balreal(sys);
[Ab,Bb,Cb,Db]=ssdata(sysb);
 
k = 2;
ak = Ab(1:k,1:k);bk=Bb(1:k,1);ck=Cb(1,1:k);
sysk = ss(ak,bk,ck,0);
sysbd = ss(0,1,0,2*sum(sv(k+1:end,1)));
z_L = impulse(sysk,t);

figure;
hold on;
plot(t,z,'b-');
plot(t,z_tilde,'r--');
plot(t,z_mod,'g');
plot(t,z_L,'k-.');
xlabel('$t$','Interpreter','latex')
ylabel('$z(t)$','Interpreter','latex')
lgd = legend('$z$','$\tilde{z}$','$z_{mod}$','$z_L$','Interpreter','latex');
lgd.FontSize = 20;
title('Impulse Responses')
hold off; 


%% Step Responses
z_step = impulse2step(z',0.01);
z_tilde_step = impulse2step(z_tilde',0.01);
z_mod_step = impulse2step(z_mod',0.01);
z_L_step = impulse2step(z_L,0.01);
figure;
hold on;
plot(t,z_step,'b-');
plot(t,z_tilde_step,'r--');
plot(t,z_mod_step,'g');
plot(t,z_L_step,'k-.');
xlabel('$t$','Interpreter','latex')
lgd = legend('$z$','$\tilde{z}$','$z_{mod}$','$z_L$','Interpreter','latex');
lgd.FontSize = 20;
title('Step Responses corresponding to the approximations')
hold off;

%% Plotting Bode amplitude diagrams
% For z
n = 1:50;
A = diag(-n.^2 * pi^2);
B = (-1).^(n)';
C = 2 * ones(1,size(n,2));
sys = ss(A,B,C,0);

[sysb,sv,T]=balreal(sys);
[Ab,Bb,Cb,Db]=ssdata(sysb);
 
k = 2;
ak = Ab(1:k,1:k);bk=Bb(1:k,1);ck=Cb(1,1:k);
sysk = ss(ak,bk,ck,0);
sysbd = ss(0,1,0,2*sum(sv(k+1:end,1)));

% For z_tilde
n = 1:8;
A = diag(-n.^2 * pi^2);
B = (-1).^(n)';
C = 2 * ones(1,size(n,2));
sys_tilde = ss(A,B,C,0);

% For z_mod
n = 1:2;
A = diag(-n.^2 * pi^2);
B = (-1).^(n)';
C = 2 * ones(1,size(n,2));
sys_mod = ss(A,B,C,0);

figure
bodeplot(sys,sys_tilde,sys_mod,sysk,{10^(-6),10^8})
lgd = legend('$h$','$\tilde{h}$','$h_{mod}$','$h_{L}$','Interpreter','latex');
lgd.FontSize = 12;
% title('Step Responses corresponding to the approximations');

%% Plotting Bode amplitude diagrams for the error systems
figure
bodeplot(sys_tilde-sys,sys_mod-sys,sysk-sys,sysbd,{10^(-6),10^8})
lgd = legend('$\tilde{h}$','$h_{mod}$','$h_{L}$','Error','Interpreter','latex');
lgd.FontSize = 12;
title('Bode diagram of error systems')
[norm(sys_tilde-sys,inf),norm(sys_mod-sys,inf),norm(sysk-sys,inf)]

function h=h(t)
    h=0;
    for k=1:50 % 50 terms considered
        h=h+(-1)^k*exp(-k^2*pi^2*t);
    end
end

function h_tilde=h_tilde(t)
    h_tilde=0;
    for k=1:8 % 8 terms considered
        h_tilde=h_tilde+(-1)^k*exp(-k^2*pi^2*t);
    end
end