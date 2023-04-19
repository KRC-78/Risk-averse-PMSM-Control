clear
close all
%%
distribution_type = "gaussian";

% ABB HDS180-4876B Parameters
f = 1e3; % Sampling frequency 1kHz
Ts = 1/f; % Sampling time 100ns
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

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*p*psi_ro);
vdss = - Lq * omegass * p*iqss;
vqss = omegass * p*psi_ro + Rs * iqss;

A = [(1 - Ts*Rs/Ld), (Ts*p*Lq*omegass/Ld), Ts*p*iqss; -(Ts*Ld*p*omegass/Lq), (1 - Ts*Rs/Lq), -Ts*p*(psi_ro/Lq + idss); 0, 1.5*Ts*p*psi_ro/Ir, (1 - Ts*Br/Ir)];
B = [Ts/Ld, 0; 0, Ts/Lq; 0, 0];

% All noise distributions are zero-mean
wbar = [0; 0; 0]; % noise mean

if distribution_type == "gaussian" % for Gaussian with stddev = 1
    tl_sd = 1; % load torque standard deviation
    M3 = [0; 0; 0];
    m4 = 1800;
    pd = makedist('Normal','mu',0,'sigma',1);
elseif distribution_type == "logistic" % for logistic with mu = 1
    tl_sd = pi / sqrt(3);
    M3 = [0; 0; 0];
    m4 = 320 * pi^4;
    pd = makedist('Logistic','mu',0,'sigma',1);
elseif distribution_type == "uniform" % for uniform from -1 to 1
    tl_sd = sqrt(1/3);
    M3 = [0; 0; 0];
    m4 = 80;
    pd = makedist('Uniform','lower',-1,'upper',1);
elseif distribution_type == "tstudent" % for t-student with nu = 5
    tl_sd = sqrt(5/3);
    M3 = [0; 0; 0];
    m4 = 900*(gamma(2.5)*gamma(0.5)*25/(sqrt(pi)*gamma(2.5)) - gamma(1.5)*gamma(1.5)*5/(sqrt(pi)*gamma(2.5)));
    pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',5);
end

W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 100000; % control sequence length

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% Risk constrained LQR
mu_nums = 10;
min_mu = -1.5; % Gaussian -3, logistic -3, uniform -2, t-student -2
max_mu = 0.5; % Gaussian 1, logistic 1, uniform 1, t-student 2
mu_array = logspace(min_mu,max_mu,mu_nums);

omega_cost = zeros(num_states, mu_nums);
omega_cost_mean = zeros(1,mu_nums);
omega_cost_stddev = zeros(1,mu_nums);
x_cost = zeros(num_states, mu_nums);
x_cost_mean = zeros(1,mu_nums);
x_cost_stddev = zeros(1,mu_nums);

for mu_count = 1:mu_nums
    mu = mu_array(mu_count);
    Q_mu = Q + 4*mu*Q*W*Q;

    V = zeros(3,3,num_states);
    K = zeros(2,3,num_states);
    S = zeros(3,3,num_states);
    T = zeros(3,3,num_states);
    l = zeros(2,num_states);
    h = zeros(2,num_states);

    V(:,:,num_states) = Q_mu;
    S(:,:,num_states) = Q;

    for n = num_states-1:-1:1
        V_prev = V(:,:,n+1);
        S_prev = S(:,:,n+1);
        T_prev = T(:,:,n+1);
        temp = inv(B'*V_prev*B+R);
        V(:,:,n) = A'*V_prev*A + Q_mu - A'*V_prev*B*temp*B'*V_prev*A;
        K(:,:,n) = -temp*B'*V_prev*A;
        S(:,:,n) = (A+B*K(:,:,n))'*S_prev+Q;
        T(:,:,n) = (A+B*K(:,:,n))'*(T_prev+V_prev);
        l(:,n) = -2*mu*temp*B'*S_prev*M3;
        h(:,n) = -temp*B'*(T_prev+V_prev)*wbar;
    end

    x = zeros(3,num_states);
    u = zeros(2,num_states);
    tl_noise = random(pd,1,num_states);

    for n = 1:num_states-1
        u(:,n) = K(:,:,n) * x(:,n) + l(:,n) + h(:,n);
        x(:,n+1) = A*x(:,n) + B*u(:,n) + [0;0;-Ts*tl_noise(n)/Ir];
        omega_cost(n, mu_count) = x(3,n)^2;
        x_cost(n,mu_count) = x(:,n)'*Q*x(:,n);
    end
    omega_cost(num_states,mu_count) = x(3, num_states)^2;
    omega_cost_mean(mu_count) = sum(omega_cost(:,mu_count))/num_states;
    omega_cost_stddev(mu_count) = sqrt(sum(omega_cost(:,mu_count).^2)/num_states - omega_cost_mean(mu_count)^2);
    x_cost(num_states,mu_count) = x(:,num_states)'*Q*x(:,num_states);
    x_cost_mean(mu_count) = sum(x_cost(:,mu_count))/num_states;
    x_cost_stddev(mu_count) = sqrt(sum(x_cost(:,mu_count).^2)/num_states - x_cost_mean(mu_count)^2);
end

%% Plot state cost mean vs omega stddev
figure
hold on
grid on
title('State cost mean vs omega cost stddev')

plot(x_cost_mean, omega_cost_stddev)

xlabel('State mean')
ylabel('Omega stddev')



