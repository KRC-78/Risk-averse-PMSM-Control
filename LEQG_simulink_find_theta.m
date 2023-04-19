clear
close all
clc
%%
% ABB HDS180-4876B Parameters
f = 1e3; % Sampling frequency 1kHz
T = 1/f; % Sampling time 100ns
omega_rpm = 1500;
omega_sync = omega_rpm*2*pi/60; % synchronous speed
Ld = 1.45e-3; % d-axis inductance 
Lq = 1.45e-3; % q-axis inductance 
Rs = 65e-3; % stator resistance
psi_ro = 0.16513; % Magnetic flux linkage 
Ir = 0.00821; % Moment of inertia
Br = 0.004; % Viscosity coefficient
p = 5; % # of pole pairs

Q = [100 0 0; 0 1 0; 0 0 30];
R = [0.8 0; 0 0.8];

tl_mean = 0; % load torque mean
tl_sd = 1; % load torque standard deviation
tl_var = tl_sd^2;

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*p*psi_ro);
vdss = - Lq * omegass * p*iqss;
vqss = omegass * p*psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*p*Lq*omegass/Ld), T*p*iqss; -(T*Ld*p*omegass/Lq), (1 - T*Rs/Lq), -T*p*(psi_ro/Lq + idss); 0, 1.5*T*p*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];
D = [0; 0; -T/Ir];

wbar = [0; 0; 0]; % noise mean
W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% LEQG
theta = -0.04; % min theta = -0.36 for stddev = 1
% Used min theta: -0.24

P = zeros(3,3,num_states);
P(:,:,num_states) = Q;
K = zeros(3,3,num_states);
K(:,:,num_states) = eye(3) - theta * Q * D * (inv(tl_var) + theta * D' * Q * D)^(-1) * D'; 
P_tilde = zeros(3,3,num_states);
P_tilde(:,:,num_states) = K(:,:,num_states) * Q;
L = zeros(2,2,num_states);
L(:,:,num_states) = -(R + B' * P_tilde(:,:,num_states) * B)^(-1);
s = zeros(3,num_states);
s(:,num_states) = -Q * xref(:,num_states);
t = zeros(1,num_states);
t(:,num_states) = xref(:,num_states)' * Q * xref(:,num_states);

for n = num_states-1:-1:1
    P_prev = P(:,:,n+1);
    P_tilde_prev = P_tilde(:,:,n+1);
    K_prev = K(:,:,n+1);
    L_prev = L(:,:,n+1);
    s_prev = s(:,n+1);
    M_prev = (inv(tl_var) + theta * D' * P_prev * D)^(-1);

    P_curr = Q + A'*P_tilde_prev*A + A'*P_tilde_prev*B*L_prev*B'*P_tilde_prev*A;
    P(:,:,n) = P_curr;
    M_curr = (inv(tl_var) + theta * D' * P_curr * D)^(-1);
    K(:,:,n) = eye(3) - theta * P_curr * D * M_curr * D';
    P_tilde(:,:,n) = K(:,:,n) * P_curr;
    L(:,:,n) = -(R + B' * P_tilde(:,:,n) * B)^(-1);
    s(:,n) = A'*K_prev*s_prev + A'*P_tilde_prev*D*tl_mean + A'*P_tilde_prev*B*L_prev*(B'*P_tilde_prev*D*tl_mean + B'*K_prev*s_prev - R*uref(:,n)) - Q*xref(:,n);
    t(n) = t(n+1) + xref(:,n)'*Q*xref(:,n) + uref(:,n)'*R*uref(:,n) + tl_mean'*D'*P_tilde_prev*D*tl_mean + 2*tl_mean'*D'*K_prev*s_prev - theta * s_prev'*D*M_prev*D'*s_prev - log(det(M_prev)/det(tl_var))/theta + (B'*P_tilde_prev*D*tl_mean + B'*K_prev*s_prev - R*uref(:,n))'*L_prev*(B'*P_tilde_prev*D*tl_mean + B'*K_prev*s_prev - R*uref(:,n));
end

x = zeros(3,num_states);
x(:,1) = xref(:,1);
u = zeros(2,num_states);
u(:,1) = uref(:,1);
tl = randn(1,num_states) * tl_sd; % stochastic noise in load torque

cost = 0;
for n = 1:num_states-1
    u(:,n) = L(:,:,n+1)*(B'*P_tilde(:,:,n+1)*(A*x(:,n) + D*tl_mean) + B'*K(:,:,n+1)*s(:,n+1) - R*uref(:,n));
    x(:,n+1) = A*x(:,n) + B*u(:,n) + D*(tl(n) - tl_mean);
    cost = cost + (x(:,n) - xref(:,n))'* Q * (x(:,n) - xref(:,n)) + (u(:,n) - uref(:,n))'* R * (u(:,n) - uref(:,n));
end

%% Plot x, u versus k
figure
hold on
grid on
title('x')

yyaxis left
plot(x(1,1:num_states-1) + xss(1))
plot(x(2,1:num_states-1) + xss(2))

xlabel('n')
ylabel('i (A)')

yyaxis right
plot(x(3,1:num_states-1) + xss(3))
ylabel('\omega (rad/s)')
legend('i_d', 'i_q','\omega')

figure
hold on
grid on
title('u')

xlabel('n')
yyaxis left
plot(u(1,1:num_states-1) + uss(1))
ylabel('v_d')

yyaxis right
plot(u(2,1:num_states-1) + uss(2))
ylabel('v_q')

legend('v_d', 'v_q')