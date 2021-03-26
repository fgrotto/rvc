clear all;
close all;

qi = 0;
qf = 10;
dot_qi = 0;
dot_qf = 0;
dot_dot_qi = 0;
dot_dot_qf = 0;
ti = 1;
tf = 8;
Ts = 0.01;
t = linspace(ti, tf, (tf-ti)/Ts);

[q, d_q, dd_q, ddd_q] = fifth_order_trajectory_siciliano(qi, dot_qi, dot_dot_qi, qf, dot_qf, dot_dot_qf, ti, tf, Ts);

plot_trajectory('Fifth polynomial trajectory (Siciliano)', t, [q; d_q; dd_q; ddd_q])

[q, d_q, dd_q, ddd_q] = fifth_order_trajectory(qi, dot_qi, dot_dot_qi, qf, dot_qf, dot_dot_qf, ti, tf, Ts);

plot_trajectory('Fifth polynomial trajectory', t, [q; d_q; dd_q; ddd_q])