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
tl_var = tl_sd^2; % load torque variance

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*psi_ro);
vdss = - Lq * omegass * iqss;
vqss = omegass * psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*Lq*omegass/Ld), T*iqss; -(T*Ld*omegass/Lq), (1 - T*Rs/Lq), -T*(psi_ro/Lq + idss); 0, 1.5*T*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];
D = [0; 0; -T/Ir];

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% LEQG controller

num_thetas = 10; % numbers of theta values
min_theta = 0;
max_theta = 0;
if tl_sd == 1
    min_theta = -1;
    max_theta = -0.1;
elseif tl_sd == 0.5
    min_theta = -3;
    max_theta = -0.3;
elseif tl_sd == 2
    min_theta = -0.22;
    max_theta = -0.02;
else 
    fprintf('Unspecified mininum theta value.')
    return;
end

theta_array = linspace(max_theta, min_theta, num_thetas);
num_trajectories = 5000; % numbers of trajectories
all_trajectory_costs = zeros(num_thetas, num_trajectories);

for theta_count = 1:num_thetas
    theta = theta_array(theta_count);
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

    for trajectory_count = 1:num_trajectories
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
        cost = cost + (x(:,num_states) - xref(:,num_states))'* Q * (x(:,num_states) - xref(:,num_states));
        all_trajectory_costs(theta_count,trajectory_count) = cost/num_states;
    end
    fprintf('Generated %d new trajectories for theta = %2.2f .\n', num_trajectories, theta);
end

if tl_sd == 0.5
    load("leqg_gaussian_05.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
elseif tl_sd == 1
    load("leqg_gaussian_1.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
elseif tl_sd == 2
    load("leqg_gaussian_2.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
end

    prev_sample_num = leqg_sample_num(end,:);
    prev_mean = leqg_mean(end,:);
    prev_stddev = leqg_stddev(end,:);
    
    curr_sample_num = zeros(1,num_thetas);
    curr_mean = zeros(1,num_thetas);
    curr_stddev = zeros(1,num_thetas);
    for theta_count = 1:num_thetas
        curr_sample_num(theta_count) = prev_sample_num(theta_count) + num_trajectories;
        curr_mean(theta_count) = (prev_mean(theta_count) * prev_sample_num(theta_count) + sum(all_trajectory_costs(theta_count,:)))/curr_sample_num(theta_count);
        curr_stddev(theta_count) = sqrt(((prev_mean(theta_count)^2 + prev_stddev(theta_count)^2)*prev_sample_num(theta_count) + sum(all_trajectory_costs(theta_count,:).^2))/curr_sample_num(theta_count) - curr_mean(theta_count)^2);
    end

    leqg_sample_num(end+1,:) = curr_sample_num;
    leqg_mean(end+1,:) = curr_mean;
    leqg_stddev(end+1,:) = curr_stddev;

    if tl_sd == 0.5
        save("leqg_gaussian_05.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
    elseif tl_sd == 1
        save("leqg_gaussian_1.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
    elseif tl_sd == 2
        save("leqg_gaussian_2.mat", "leqg_stddev", "leqg_mean", "leqg_sample_num");
    end

 plot_lqr_leqg