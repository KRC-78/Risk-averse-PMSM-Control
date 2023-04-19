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

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*p*psi_ro);
vdss = - Lq * omegass * p*iqss;
vqss = omegass * p*psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*p*Lq*omegass/Ld), T*p*iqss; -(T*Ld*p*omegass/Lq), (1 - T*Rs/Lq), -T*p*(psi_ro/Lq + idss); 0, 1.5*T*p*psi_ro/Ir, (1 - T*Br/Ir)];
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


controller_gain = L(:,:,1)*B'*P(:,:,1)*A;
controller_offset = L(:,:,1)*B'*s(:,1);