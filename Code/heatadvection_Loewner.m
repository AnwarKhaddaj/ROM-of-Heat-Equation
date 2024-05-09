% Close all previous figures and clear all variables
close all;
clear all;

% Discretized Nth order system
N_disc = 500;

eta = 0.01;
h_disc = 1/(N_disc+1);

A = spdiags([1+2*h_disc*eta*ones(N_disc, 1), -2*(1+h_disc*eta)*ones(N_disc, 1), ones(N_disc, 1)], -1:1, N_disc, N_disc);
A = full(A);

A(1,1) = -1;
A(end, end) = -1-2*h_disc*eta;
A = (1/h_disc^2)*A;

B = zeros(N_disc, 1);
B(end, 1) = 1;
B = 1/h_disc * B;

C = zeros(1, N_disc);
C(1, 1) = 1;

E = eye(N_disc, N_disc);
sys = ss(A,B,C,0);

% figure;
% % Plot impulse response
% hold on;
% t = linspace(0,1,1000);
% 
% plot(T,IRorig,'b-');
% % plot(T,IRk,'g-.'); 
% legend('original','Loewner')
% title('Impulse response of reduced system with Loewner')
% hold off;
% 
% keyboard

%%%%% Reduced Order Modeling : Loewner
N = 250; % number of interpolation points
% om = linspace(1,250,N);
om = logspace(-1,2,N); %(-1,2)
lam = kron(om(1,1:2:N),[1j -1j]);
mu = kron(om(1,2:2:N), [1j -1j]).';

R = ones(1,N); L=R';
W = zeros(1,size(mu,1));
V = zeros(size(mu,1),1);
for i=1:size(mu,1)
    W(i) = Z(C,E,A,B,lam(i));
    V(i) = Z(C,E,A,B,mu(i));
end

% left anf right interpolation data

L_c =(V*R-L*W)./(mu*R-L*lam);
SL_c=((mu.*V)*R-L*(W.*lam))./(mu*R-L*lam);

% Converting Loewner and Shifted Loewner matrix to real form
J0=[1 -1i;1 1i]; J=kron(eye(N/2),J0);
L_r = J'*L_c*J;
SL_r = J'*SL_c*J;
Vr = J'*V; Lr = J'*L;
Wr = W*J; Rr = R*J;

% SVD of Loewner pencil and its transpose and its resulting projection
% matrices x and y
[u1,s1,v1] = svd([L_r SL_r]);
[u2,s2,v2] = svd([L_r; SL_r]);
figure;
loglog(diag(s1))
title("Singular Values of Loewner Pencil")
k = 183; %37
x=u1(:,1:k);
y=v2(:,1:k);

% e1=x'*L_r*y;
% a1=x'*SL_r*y;
% b1=x'*Vr;
% c1=Wr*y;
% sysL=ss(e1,a1,b1,c1);
Qk = x'*SL_r*y; Pk = x'*L_r*y; Wk = Wr*y; Vk = x'*Vr; Lk=x'*Lr; Rk=Rr*y;
sysLk=dss(-Qk,Vk,Wk,0,-Pk); 

% e1=x'*L_r*y;a1=x'*SL_r*y;b1=x'*Vr;c1=Wr*y;
% sysLk=ss(e1\a1,e1\b1,c1,0);

figure;
sigma(sysLk)
t = linspace(0,1,1000); % time step = 0.1 on [0,1]

figure;
% Plot impulse response
hold on;
s = linspace(0,1,1000);
[IRorig,T] = impulse(sys);
plot(T,IRorig,'b-');
[IRk,T]=impulse(sysLk,T);
plot(T,IRk,'g-.'); 
legend('original','Loewner')
title('Impulse response of reduced system with Loewner')
hold off;

% n = 0:50;
% A = diag(-n.^2 * pi^2);
% B = (-1).^(n)';
% C = [1, 2 * ones(1,size(n,2)-1)];

figure;
poles_orig = pole(sys);
poles = pole(sysLk);
hold on;
% plot(real(poles_orig),imag(poles_orig),".")
plot(real(poles),imag(poles),"o")
title('Poles of the reduced system with Loewner')
xlabel("Re(z)")
ylabel("Im(z)")
hold off;

figure
bodeplot(sys, {10^(-8),10^(2)}, 'b')
hold on;
bodeplot(sysLk,{10^(-8),10^(2)}, 'r-.')
bodeplot(sys-sysLk,{10^(-8),10^(2)}, 'g-.')
lgd = legend('$z$','$z_{Loewner}$','Error','Interpreter','latex');
lgd.FontSize = 12;
hold off;

% function Z=Z(C,E,A,B,s)
%     eta = 0.01;
%     Z = (1/(s.*exp(eta))).*(sqrt(eta^2+s)/sinh(sqrt(eta^2+s)));
% %     Z = 1./((s*e^)sqrt(s).*sinh(sqrt(s)));
% end

function Z=Z(C,E,A,B,t)
    Z = C*((t*E-A)\B);
end

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