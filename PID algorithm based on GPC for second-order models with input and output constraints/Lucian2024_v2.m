% Paper simulation
% https://doi.org/10.1016/j.ifacol.2024.08.002

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% 2nd order model
% Continuous time
s = tf('s');
Ps = (-0.8*s+1)/(1.5*s+1)^2

% Sampling time
Ts = 0.1;

% Discrete time
Pz = c2d(Ps,Ts,'zoh')
[num,den] = tfdata(Pz,'v');
% num = round(num, 3, 'decimals');
% den = round(den, 3, 'decimals');
% Pz = tf(num,den,Ts)

%% Controller tuning parameters
% Prediction horizon
N = 20;
% Control horizon
Nu = 1;

% Weight values
lambda = 0;
lambda_e = 1000;
Q = eye(N);
for aux=1:8
    Q(aux,aux) = 0;
end

% Constraints
du_min = -0.5;
du_max = 0.5;
u_min = 0;
u_max = 0.9;
y_min = 0;
y_max = 0.7;

%% CARIMA model
Delta = tf([1 -1],[1],Ts,'Variable','z^-1'); % Δ = 1 - z^-1
[delta, ~] = tfdata(Delta, 'v');

% Polynomials b(z^-1) and a(z^-1)
[b,a] = tfdata(Pz,'v');

% Disturbance rejection polynomial c(z^-1) 
alf = 0.5;
c  = conv([1 -alf],[1 -alf]);

a_til = conv(delta,a); % a_til = Δ*a

Max = max([length(a_til)-1,length(b)-1,length(c)-1]);
a_til = [a_til zeros(1,Max-(length(a_til)-1))]
b = [b zeros(1,Max-(length(b)-1))]
c = [c zeros(1,Max-(length(c)-1))]

% State-Space cannonical observable form
A = [-a_til(2:end)' [eye(length(a_til)-2); zeros(1,length(a_til)-2)]]
B = [b(2:end)']
D = [c(2:end)'-a_til(2:end)']
H = [1 zeros(1,Max-1)]

% Prediction Matrix
for i=1:N
   G(i,:) = H*A^(i-1)*B;
end
for j=2:N
    G = [G [0; G(1:end-1,j-1)]];
end
% % for k=Nu+1:N
% %     G(:,Nu) = G(:,Nu) + G(:,k);
% % end
G = G(:,1:Nu)

for i=1:N
    F(i,:) = H*A^i;
end
F

for i=1:N
    E(i,1) = H*A^(i-1)*D;
end
E

%% Simulation parameters
% Time parameters
Tsim = 60;
qntd_k = ceil(Tsim/Ts);

% Reference parameters
r1 = 0.5;
Tr1 = ceil(1/Ts);
r2 = 0.3;
Tr2 = ceil(45/Ts);

% Disturbance parameters
pert = -0.1;
Tpert = ceil(30/Ts);

% Reference and disturbance arrays
ref = [zeros(Tr1,1); r1*ones(Tr2-Tr1,1); (r1+r2)*ones(qntd_k-Tr2,1)];
q = [zeros(Tpert,1); pert*ones(qntd_k-Tpert,1)];

%return
%% PID without constraints
% unconstrained controller gains
K = inv(G'*Q*G+lambda)*G'*Q;
K = K(1,:);
Kr = sum(K);
KF = K*F;
KE = K*E;
% KF = acker(A,B,[1 1 1]*0.7)
% Kr = inv(H*inv(eye(3)-A+B*KF)*B)
% KE = acker(A,B,eig((A)-B*KF)'); KE = KE(1);

% augmented system
Ac = A-D*H+B*(KE*H-KF);
Bc = [B*Kr,(D-B*KE)];
Cc = KE*H-KF;
Dc = [Kr,-KE];

Controlador = tf(ss(Ac,Bc,Cc,Dc,Ts));
[nC1,dC1] = tfdata(Controlador(1),'v');
[nC2,dC2] = tfdata(Controlador(2),'v');

C_PID =  -zpk(minreal(tf(nC2,conv(dC2(1:end-1),delta),Ts,'Variable','z^-1')))
F_PID =  -zpk(minreal(tf([nC1 0],nC2,Ts,'Variable','z^-1')))

%% Simulation
tic
fprintf('Simulation Runtime: ')
out = sim('simu_Lucian2024.slx');
t = out.t(2:end);
y_GPC = out.y_GPC(2:end);
u_GPC = out.u_GPC(2:end);
du_GPC = out.du_GPC(2:end);
y_PID = out.y_PID(2:end);
u_PID = out.u_PID(2:end);
du_PID = out.du_PID(2:end);
y_PIDmapping = out.y_PIDmapping(2:end);
u_PIDmapping = out.u_PIDmapping(2:end);
du_PIDmapping = out.du_PIDmapping(2:end);
toc

figure
subplot(3,1,1)
hold on
stairs(t,y_GPC, 'm', 'linewidth', 2)
stairs(t,y_PIDmapping, 'b', 'linewidth', 2)
stairs(t,y_PID, '-.g', 'linewidth', 2)
stairs(t,ref, '--r', 'linewidth', 2)
stairs(t,y_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,y_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('y')
legend('GPC', 'PID with mapping', 'Proposed PID', 'Reference', 'Constraints', 'location', 'best')
title('Output')
subplot(3,1,2)
hold on
stairs(t,u_GPC, 'm', 'linewidth', 2)
stairs(t,u_PIDmapping, 'b', 'linewidth', 2)
stairs(t,u_PID, '-.g', 'linewidth', 2)
stairs(t,u_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,u_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('u')
title('Control signal')
subplot(3,1,3)
hold on
stairs(t,du_GPC, 'm', 'linewidth', 2)
stairs(t,du_PIDmapping, 'b', 'linewidth', 2)
stairs(t,du_PID, '-.g', 'linewidth', 2)
stairs(t,du_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,du_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('$\Delta u$')
title('Control signal rate')
sgtitle("\textbf{Output and Control Signals for N = }" + N + "\textbf{ and Nu = }" + Nu)

fprintf('\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')

% fim