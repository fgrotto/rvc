clear all;
close all;

addpath('../ur5_model');

q_s = [0 2 -2];
t_s = [0 2 4];
maxTau = [9 13 10 9 9 3];
Ts = 0.001;

dh = load_dh_parameters;

[q,dq,ddq,tau,tau_raw, lambda,q_raw,dq_raw,ddq_raw] = model_based_trajectory(dh,q_s,t_s,maxTau, 'acc', Ts);

for i = 1:length(maxTau)
    disp(['q',num2str(i),' tau max obtained after rescaling ', num2str(max(abs(tau(i,:))))]);
end

plot_model_based_trajectory(q,dq,ddq,Ts,q_s,t_s./lambda)
plot_torques(tau,t_s./lambda, Ts, maxTau, "new torques")
plot_torques(tau_raw,t_s./lambda, Ts, maxTau, "old torques")
