%% Multipoint trajectory with discontinuos acceleration
close all; 
clear all;

q_d = [10 20 0 30 40];
t_d = [0 2 4 8 10];
d_qi = -3;
d_qf = 3;
Ts = 0.01;

[q,d_q,dd_q] = multipoint_trajectory(q_d, t_d, d_qi, d_qf, Ts);
plot_trajectory('Multipoint trajectory', t_d(1):Ts:t_d(end), [q; d_q; dd_q])

%% Multipoint trajectory with continuous acceleration
close all; 
clear all;

q_d = [10 20 0 30 40];
t_d = [0 2 4 8 10];
d_qi = 0;
d_qf = 0;
Ts = 0.01;

[q,d_q,dd_q] = multipoint_trajectory_acc_cont(q_d, t_d, d_qi, d_qf, Ts);
plot_trajectory('Multipoint trajectory', t_d(1):Ts:t_d(end), [q; d_q; dd_q])
