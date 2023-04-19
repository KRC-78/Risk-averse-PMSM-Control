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
tl_sd = 2; % load torque standard deviation

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*psi_ro);
vdss = - Lq * omegass * iqss;
vqss = omegass * psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*Lq*omegass/Ld), T*iqss; -(T*Ld*omegass/Lq), (1 - T*Rs/Lq), -T*(psi_ro/Lq + idss); 0, 1.5*T*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];

W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% LQG controller

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
    s(:,n) = A'*(s_prev + P_prev * B*L_prev*(B'*s_prev - R*uref(:,n))) - Q*xref(:,n);
    t(n) = t(n+1) + xref(:,n)'*Q*xref(:,n) + uref(:,n)'*R*uref(:,n) + trace(P_prev*W) + (B'*s_prev - R*uref(:,n))'*L_prev*(B'*s_prev - R*uref(:,n));
end

%% M trajectories

for count = 1:19
num_trajectories = 10000;
all_trajectory_costs = zeros(1,num_trajectories);

for trajectory_count = 1:num_trajectories
    x = zeros(3,num_states);
    x(:,1) = xref(:,1);
    u = zeros(2,num_states);
    u(:,1) = uref(:,1);
    tl = randn(1,num_states) * tl_sd; % stochastic noise in load torque

    cost = 0;
    for n = 1:num_states-1
        u(:,n) = L(:,:,n+1)*(B'*P(:,:,n+1)*A*x(:,n) + B'*s(:,n+1) - R*uref(:,n));
        x(:,n+1) = A*x(:,n) + B*u(:,n) + [0;0;-T*(tl(n)-tl_mean)/Ir];
        cost = cost + (x(:,n) - xref(:,n))'* Q * (x(:,n) - xref(:,n)) + (u(:,n) - uref(:,n))'* R * (u(:,n) - uref(:,n));
    end
    cost = cost + (x(:,num_states) - xref(:,num_states))'* Q * (x(:,num_states) - xref(:,num_states));
    all_trajectory_costs(trajectory_count) = cost/num_states;
end

if tl_sd == 0.5
    load("lqr_gaussian_05.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
elseif tl_sd == 1
    load("lqr_gaussian_1.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
elseif tl_sd == 2
    load("lqr_gaussian_2.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
end

prev_sample_num = lqr_sample_num(end);
prev_mean = lqr_mean(end);
prev_stddev = lqr_stddev(end);

curr_sample_num = prev_sample_num + num_trajectories;
curr_mean = (prev_mean * prev_sample_num + sum(all_trajectory_costs))/curr_sample_num;
curr_stddev = sqrt(((prev_mean^2 + prev_stddev^2)*prev_sample_num + sum(all_trajectory_costs.^2))/curr_sample_num - curr_mean^2);

lqr_sample_num(end+1) = curr_sample_num;
lqr_mean(end+1) = curr_mean;
lqr_stddev(end+1) = curr_stddev;

if tl_sd == 0.5
    save("lqr_gaussian_05.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
elseif tl_sd == 1
    save("lqr_gaussian_1.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
elseif tl_sd == 2
    save("lqr_gaussian_2.mat", "lqr_stddev", "lqr_mean", "lqr_sample_num");
end
fprintf('%d total trajectories.', curr_sample_num)

end

fprintf('Generated %d new trajectories. Total: %d, mean: %4.3f, stddev: %4.3f .\n', num_trajectories, curr_sample_num, curr_mean, curr_stddev);