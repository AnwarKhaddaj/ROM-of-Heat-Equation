%  Reduced Models of Heat Diffusion Equation
%  ELEC519 Final Project
%
%  AUTHORS:  Karthik Goli, Anwar Khaddaj, and Carlos Taveras
%            Rice University
%            April 12, 2024
%            Last Modified by Anwar: April 16, 2024 3:20pm

% Close all previous figures and clear all variables
close all; clear;

A_func = @(mm, etaa) diag(-mm.^2 * pi^2 -etaa^2);
B_func = @(mm) (-1).^mm';
C_func = @(mm, etaa) (1 + etaa^2./(mm.^2 * pi^2));
sys_func = @(mm, etaa) ss(A_func(mm, etaa), B_func(mm), C_func(mm, etaa), 0);
h_func   = @(tt, mm, etaa) impulse(tt, sys_func(mm, etaa));
z_func = @(tt, mm, etaa) (etaa/sinh(etaa) + (2 * h_func(tt, mm, etaa) ) ) * exp(-etaa);

%% Plotting impulse responses z(t)
t = linspace(0,1,1000); % time step = 0.1 on [0,1]

n = 1:50; eta = 0.01;
n_mod = 1:2; n_tild = 1:8;

sys = sys_func(1:50, eta);
sys_tilde = sys_func(1:8, eta);
sys_mod = sys_func(1:2, eta);

[sysb,sv,T]=balreal(sys); [Ab,Bb,Cb,Db]=ssdata(sysb);
k = 2; ak = Ab(1:k,1:k);bk=Bb(1:k,1);ck=Cb(1,1:k);
sys_bt = ss(ak,bk,ck,0); sysbd = ss(0,1,0,2*sum(sv(k+1:end,1)));


z = z_func(t, 1:50, eta);
z_mod = z_func(t, 1:2, eta);
z_tilde = z_func(t, 1:8, eta);
z_L = (eta/sinh(eta) + 2 * impulse(sys_bt,t)) *  exp(-eta);

% Plotting poles 
poles = pole(sys);
figure;
plot(real(poles),imag(poles),'.');
title('Poles of the original system');

% impulse response
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
z_step = impulse2step(z,0.01);
z_tilde_step = impulse2step(z_tilde,0.01);
z_mod_step = impulse2step(z_mod,0.01);
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

%% Bode Plots 

figure;
bodeplot(sys,sys_tilde,sys_mod,sys_bt,{10^(-6),10^8})
lgd = legend('$z$','$\tilde{z}$','$z_{mod}$','$z_{L}$','Interpreter','latex');
lgd.FontSize = 12;
% title('Step Responses corresponding to the approximations');

%% Plotting Bode amplitude diagrams for the error systems
figure
bodeplot(sys_tilde-sys,sys_mod-sys,sys_bt-sys,sysbd,{10^(-6),10^8})
lgd = legend('$\tilde{z}$','$z_{mod}$','$z_{L}$','Error','Interpreter','latex');
lgd.FontSize = 12;
title('Bode diagram of error systems')

