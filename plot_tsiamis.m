distribution_type = "uniform";

if distribution_type == "gaussian" % for Gaussian with stddev = 1
    load('tsiamis_gaussian.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
elseif distribution_type == "logistic" % for logistic with mu = 1
    load('tsiamis_logistic.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
elseif distribution_type == "uniform" % for uniform from -1 to 1
    load('tsiamis_uniform.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
elseif distribution_type == "tstudent" % for t-student with nu = 5
    load('tsiamis_tstudent.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
end

%% Plot state cost mean vs omega stddev
figure ('position',[100 100 1600 300])
subplot(1,3,3)
plot(x_cost_mean, omega_cost_stddev, 'bo--')
grid on

subplot(1,3,2)
semilogx(mu_array, x_cost_mean, 'bo--')
grid on

subplot(1,3,1)
semilogx(mu_array, omega_cost_stddev, 'bo--')
grid on