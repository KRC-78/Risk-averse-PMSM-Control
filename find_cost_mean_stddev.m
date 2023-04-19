% clear all
% 
% load('data.mat')
% 
% theta = -0.45;
% 
% id = ans(2,:);
% iq = ans(3,:);
% omega = ans(4,:);
% vd = ans(5,:);
% vq = ans(6,:);
% 
% save('leqg_045_simulink.mat')
% %% 
% len = length(id);
% 
% cost = zeros(1,len);
% x = [id; iq; omega];
% u = [vd; vq];
% Q = diag([100, 1, 30]);
% R = diag([0.8,0.8]);
% for count = 1:len
%     cost(count) = x(:,count)'*Q*x(:,count) + u(:,count)'*R*u(:,count);
% end
% 
% curr_mean = sum(cost)/len
% curr_stddev = sqrt(sum(cost.^2)/len - curr_mean^2)
% %%
% load('simulink_leqg_lqr.mat')
% 
% theta_array(end+1) = theta;
% leqg_mean(end+1) = curr_mean;
% leqg_stddev(end+1) = curr_stddev;
% 
% save('simulink_leqg_lqr.mat','theta_array', 'leqg_stddev', 'leqg_mean')

%% 
load('simulink_leqg_lqr.mat')
figure ('position',[100 100 1600 300])
subplot(1,3,1)
plot(theta_array(1,2:end), leqg_mean(1,2:end), 'bo--')
hold on
plot(theta_array(1,1), leqg_mean(1,1), 'r*')
grid on
xlim([-0.5 0.05])

subplot(1,3,2)
plot(theta_array(1,2:end), leqg_stddev(1,2:end), 'bo--')
hold on
plot(theta_array(1,1), leqg_stddev(1,1), 'r*')
grid on
xlim([-0.5 0.05])

subplot(1,3,3)
plot(leqg_mean(1,2:end), leqg_stddev(1,2:end), 'bo--')
hold on
plot(leqg_mean(1,1), leqg_stddev(1,1), 'r*')
grid on