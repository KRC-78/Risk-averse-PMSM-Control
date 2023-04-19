clear
close all
clc
%%
% Parameters
f = 10e3; % Sampling frequency 10kHz
T = 1/f; % Sampling time 100ns
omega_sync = 314.159; % synchronous speed 314.159rad/s (3000 rpm)
Ld = 0.8e-3; % d-axis inductance 0.8mH
Lq = 0.8e-3; % q-axis inductance 0.8mH
Rs = 45e-3; % stator resistance 45mΩ
psi_ro = 0.65; % Magnetic flux linkage 0.65Wb
Ir = 0.35; % Moment of inertia 0.35 kg.m^2
Br = 0.004; % Viscosity coefficient

Q = [100 0 0; 0 1 0; 0 0 30];
R = [0.8 0; 0 0.8];

tl_mean = 0; % load torque mean
tl_sd = 1; % load torque standard deviation

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*psi_ro);
vdss = - Lq * omegass * iqss;
vqss = omegass * psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*Lq*omegass/Ld), T*iqss; -(T*Ld*omegass/Lq), (1 - T*Rs/Lq), -T*(psi_ro/Lq + idss); 0, 1.5*T*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];

wbar = [0; 0; 0]; % noise mean
W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% LQG

P = zeros(3,3,num_states);
P(:,:,num_states) = Q;
s = zeros(3,num_states);
s(:,num_states) = -Q*xref(:,num_states);
t = zeros(1,num_states);
t(:,num_states) = xref(:,num_states)'*Q*xref(:,num_states);
L = zeros(2,2,num_states);
L(:,:,num_states) = -(R+B'*Q*B)^(-1);

for n = num_states-1:-1:1
    P_prev = P(:,:,n+1);
    s_prev = s(:,n+1);
    L_prev = L(:,:,n+1);
    P(:,:,n) = Q + A'*P_prev*A + A'*P_prev*B*L_prev*B'*P_prev*A;
    L(:,:,n) = -(R+B'*P(:,:,n)*B)^(-1);
    s(:,n) = A'*(s_prev + P_prev * (wbar + B*L_prev*(B'*P_prev*wbar + B'*s_prev - R*uref(:,n)))) - Q*xref(:,n);
    t(n) = t(n+1) + xref(:,n)'*Q*xref(:,n) + uref(:,n)'*R*uref(:,n) + 2*s_prev'*wbar + trace(P_prev*W) + (B'*P_prev*wbar + B'*s_prev - R*uref(:,n))'*L_prev*(B'*P_prev*wbar + B'*s_prev - R*uref(:,n));
end

x = zeros(3,num_states);
x(:,1) = xref(:,1);
u = zeros(2,num_states);
u(:,1) = uref(:,1);
tl_noise = randn(1,num_states) * tl_sd; % stochastic noise in load torque

cost = 0;
for n = 1:num_states-1
    u(:,n) = L(:,:,n+1)*(B'*P(:,:,n+1)*(A*x(:,n)+wbar) + B'*s(:,n+1) - R*uref(:,n));
    x(:,n+1) = A*x(:,n) + B*u(:,n) + wbar + [0;0;-T*tl_noise(n)/Ir];
    cost = cost + (x(:,n) - xref(:,n))'* Q * (x(:,n) - xref(:,n)) + (u(:,n) - uref(:,n))'* R * (u(:,n) - uref(:,n));
end
cost = cost + (x(:,num_states) - xref(:,num_states))'* Q * (x(:,num_states) - xref(:,num_states));
cost/num_states
%% Plot x, u versus k
figure
hold on
grid on
title('x')

yyaxis left
plot(x(1,1:num_states-1)+xss(1))
plot(x(2,1:num_states-1)+xss(2))

xlabel('n')
ylabel('i (A)')

yyaxis right
plot(x(3,1:num_states-1)+xss(3))
ylabel('\omega (rad/s)')
legend('i_d', 'i_q','\omega')

figure
hold on
grid on
title('u')

xlabel('n')
yyaxis left
plot(u(1,1:num_states-1)+uss(1))
ylabel('v_d')

yyaxis right
plot(u(2,1:num_states-1)+uss(2))
ylabel('v_q')

legend('v_d', 'v_q')