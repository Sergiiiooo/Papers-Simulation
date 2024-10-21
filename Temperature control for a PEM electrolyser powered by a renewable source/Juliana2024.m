% Paper simulation
% https://doi.org/10.17979/ja-cea.2024.45.10894

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% Simulation parameters
% Prediction horizon
N = 10;
% Control horizon
Nu = 10;

% Weight gains
Q = 1*eye(N);
R = 1*eye(Nu);
S = 0*eye(Nu);

% Sampling time
Ts = 100;

% Temperatura ambiente
cond_inicial = 0;

% Simulation parameters
Tsim = 3600*24;
Tsim = 100000;
Tref = 0;
ref = 1; ref = ref-cond_inicial;

%% Model
% Continuos parameters
s = tf('s');
ac = 9.6832*10^-4;
bc = -1.0482*10^-4;
m1c = 2.3090*10^-4;
m2c = 9.5293*10^-4;

% x(k+1) = ad x(k) + bd u(k) + m1d d1(k) + m2d d2(k)
ad = 0.9077; a = [1 -ad];
bd = -0.01; b = [0 bd];
m1d = 0.022;
m2d = 0.0908;

Gc = tf(b,[1 -ad],Ts,'Variable','z^-1');
G2 = tf(m1d,[1 -ad],Ts,'Variable','z^-1');
G3 = tf(m2d,[1 -ad],Ts,'Variable','z^-1');

% Delta
Delta = tf([1 -1],[1],Ts,'Variable','z^-1'); % Δ = 1 - z^-1
[delta, ~] = tfdata(Delta, 'v');

% Disturbance rejection polynomial c(z^-1) 
c  = [1 -1.62 0.6561];
c = conv([1 -0.5],[1 -0.5]);

a_til = conv(delta,a); % a_til = Δ*a

Max = max([length(a_til)-1,length(b)-1,length(c)-1]);
a_til = [a_til zeros(1,Max-(length(a_til)-1))];
b = [b zeros(1,Max-(length(b)-1))];
c = [c zeros(1,Max-(length(c)-1))];

% State-Space cannonical observable form
A = [-a_til(2:end)' [eye(length(a_til)-2); zeros(1,length(a_til)-2)]]
D = [c(2:end)'-a_til(2:end)']
C = [1 zeros(1,Max-1)]

% Prediction Matrix
for i=1:N
   G(i,:) = ad^(i-1)*bd;
end
for j=2:N
    G = [G [0; G(1:end-1,j-1)]];
end
% % for k=Nu+1:N
% %     G(:,Nu) = G(:,Nu) + G(:,k);
% % end
G = G(:,1:Nu)

for i=1:N
    F1(i,:) = ad^i;
    F2(i,:) = F1(i,:);
    F3(i,:) = F1(i,:);
    F4(i,:) = C*A^i;
end
F1
F2
F3
F4

for i=1:N
    E(i,1) = C*A^(i-1)*D;
end
E

H1 = zeros(N,1); H2 = zeros(N,1);
for i = 1:N
    for j = 0:i-1
        H1(i) = H1(i) + ad^j * m1d;
        H2(i) = H2(i) + ad^j * m2d;
    end
end
H1
H2

%% GPC controller tuning
M = inv(tril(ones(N)));
T = eye(N,1);

K = inv(G'*Q*G + M'*R*M + S)*G'*Q;
K = K(1,:)

kr = sum(K)

ke = inv(G'*Q*G + M'*R*M + S)*M'*R*T;
ke = ke(1,:)

K1 = K*H1
K2 = K*H2
K3 = K*F1
K4 = K*E - ke
K5 = K*(F4 - E*C) + ke*C

z = tf('z',Ts);
V = K5*inv(z*eye(length(A)) - A + D*C)*D + K4
% return
%% Simulation and Figures
out = sim('simu_Juliana2024.slx');
y = out.y+cond_inicial;
r = out.r+cond_inicial;
u = out.u;
t = out.t; t = t/3600;
d1 = out.d1;
d2 = out.d2;

figure
subplot(2,1,1)
hold on
stairs(t,y,'r','LineWidth',2)
stairs(t,r,'--k','LineWidth',2)
legend('$$T_{el}$$', 'Ref', 'location', 'best')
grid on
axis tight
xlabel('Time (h)')
ylabel('Temperature (ºC)')
title('Electrolyser Temperature')
subplot(2,1,2)
hold on
stairs(t,u,'b','LineWidth',2)
legend('$$\dot{Q}_{cool}$$', 'location', 'best')
grid on
axis tight
xlabel('Time (h)')
ylabel('Heat (W)')
title('Heat dissipated with the Cooling System')

figure
subplot(2,1,1)
hold on
stairs(t,d1,'g','LineWidth',2)
legend('$$I_{el}$$', 'location', 'best')
grid on
axis tight
xlabel('Time (h)')
ylabel('Current (A)')
title('Electrolyser Current supplied by PV')
subplot(2,1,2)
hold on
stairs(t,d2,'m','LineWidth',2)
legend('$$T_{amb}$$', 'location', 'best')
grid on
axis tight
xlabel('Time (h)')
ylabel('Temperature (ºC)')
title('Ambient Temperature')




fprintf('\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')
% fim