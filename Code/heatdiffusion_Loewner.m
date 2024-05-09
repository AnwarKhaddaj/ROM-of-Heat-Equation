% Close all previous figures and clear all variables
close all;
clear all;

%%%%% Reduced Order Modeling : Loewner
N = 250; % number of interpolation points
% om = linspace(1,250,N);
om = logspace(-1,2,N); %(-1,2)
lam = kron(om(1,1:2:N),[1j -1j]);
mu = kron(om(1,2:2:N), [1j -1j]).';

R = ones(1,N); L=R';
W = Z(lam) ; V = Z(mu); % left anf right interpolation data

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
k = 37;
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
t = linspace(0,1,1000); 
[IRk,T]=impulse(sysLk,t);
figure;
% Plot impulse response
hold on;
t = linspace(0,1,1000); % time step = 0.1 on [0,1]
z = 1+2*h(t);
plot(t,z,'b-');
plot(T,IRk,'g-.'); 
legend('original','Loewner')
title('Impulse response of reduced system with Loewner')
hold off;

n = 0:50;
A = diag(-n.^2 * pi^2);
B = (-1).^(n)';
C = [1, 2 * ones(1,size(n,2)-1)];

sys = ss(A,B,C,0);

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
bodeplot(sys,sysLk)
hold on;
bodeplot(sys-sysLk)
lgd = legend('$sys$','$sys_{Loewner}$','Error','Interpreter','latex');
lgd.FontSize = 12;
hold off;

function Z=Z(s)
    Z = 1./(sqrt(s).*sinh(sqrt(s)));
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